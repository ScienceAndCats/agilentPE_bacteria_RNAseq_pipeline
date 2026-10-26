"""Bounded, streaming I/O for native featureCounts tab-separated tables.

This is literal TSV, not quoted CSV: quotes and semicolon-separated annotation
lists are preserved. The 8 MiB limit applies to each complete UTF-8 physical
line, including its newline, before decoding or splitting. Files are never
loaded in full. Counts and feature lengths must fit an unsigned 64-bit integer.
"""

from contextlib import contextmanager
from pathlib import Path
import re


ANNOTATION_COLUMNS = ("Geneid", "Chr", "Start", "End", "Strand", "Length")
MAX_LINE_BYTES = 8 * 1024 * 1024
MAX_COUNT = (1 << 64) - 1
MAX_SUMMARY_ROWS = 256
MAX_SAMPLE_COLUMNS = 4096
CONTROL = re.compile(rb"[\x00-\x08\x0b-\x1f\x7f]")


def bounded_lines(handle, source, *, max_line_bytes=MAX_LINE_BYTES):
    """Read at most limit+1 bytes, even for a missing newline or bad comment."""
    if type(max_line_bytes) is not int or max_line_bytes < 1:
        raise ValueError("featureCounts line limit must be a positive integer")
    number = 0
    while True:
        raw = handle.readline(max_line_bytes + 1)
        if not raw:
            return
        number += 1
        if len(raw) > max_line_bytes:
            raise ValueError(f"featureCounts line exceeds {max_line_bytes} bytes at {source}:{number}")
        raw = raw.removesuffix(b"\n").removesuffix(b"\r")
        if CONTROL.search(raw):
            raise ValueError(f"Invalid control character in featureCounts table at {source}:{number}")
        try:
            yield number, raw.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise ValueError(f"Invalid UTF-8 in featureCounts table at {source}:{number}") from exc


def count_value(value, source, number, label="count"):
    # Check length before int() so malformed million-digit numbers cannot cause
    # an expensive conversion or an exception from Python's own integer limit.
    if (re.fullmatch(r"[0-9]{1,20}", value) is None or int(value) > MAX_COUNT):
        raise ValueError(f"Invalid featureCounts {label} at {source}:{number}: {value[:40]!r}")
    return value


def _table_rows(handle, source, annotation_columns, max_line_bytes):
    header = None
    summary = annotation_columns == ("Status",)
    seen_statuses = set()
    for number, line in bounded_lines(handle, source, max_line_bytes=max_line_bytes):
        if not line or line.startswith("#"):
            continue
        width = len(header) if header is not None else len(annotation_columns) + MAX_SAMPLE_COLUMNS
        row = line.split("\t", width)
        if header is None:
            if (tuple(row[:len(annotation_columns)]) != annotation_columns
                    or not len(annotation_columns) < len(row) <= width):
                raise ValueError(f"Unexpected featureCounts header at {source}:{number}")
            samples = row[len(annotation_columns):]
            if any(not name for name in samples) or len(samples) != len(set(samples)):
                raise ValueError(f"Duplicate or empty featureCounts sample columns at {source}:{number}")
            header = row
            yield row
            continue
        if len(row) != len(header):
            raise ValueError(f"Unexpected featureCounts column count at {source}:{number}: "
                             f"expected {len(header)}, found {len(row)}")
        if any(not field for field in row[:len(annotation_columns)]):
            raise ValueError(f"Empty featureCounts annotation at {source}:{number}")
        if summary:
            if re.fullmatch(r"[A-Za-z][A-Za-z0-9_]{0,255}", row[0]) is None:
                raise ValueError(f"Invalid featureCounts status at {source}:{number}")
            if row[0] in seen_statuses:
                raise ValueError(f"Mismatched or duplicate featureCounts status at {source}:{number}")
            if len(seen_statuses) >= MAX_SUMMARY_ROWS:
                raise ValueError(f"Too many featureCounts summary rows at {source}:{number}")
            seen_statuses.add(row[0])
        else:
            count_value(row[5], source, number, "feature length")
        for value in row[len(annotation_columns):]:
            count_value(value, source, number)
        yield row
    if header is None:
        raise ValueError(f"Missing featureCounts table header: {source}")
    if summary and "Assigned" not in seen_statuses:
        raise ValueError(f"No Assigned status in featureCounts summary: {source}")


@contextmanager
def open_table(path, annotation_columns=ANNOTATION_COLUMNS, *, max_line_bytes=MAX_LINE_BYTES):
    """Yield validated header and data rows; close the file even on failure."""
    path = Path(path)
    columns = tuple(annotation_columns)
    if columns not in (ANNOTATION_COLUMNS, ("Status",)):
        raise ValueError("Expected a featureCounts matrix or assignment summary")
    with path.open("rb") as handle:
        yield _table_rows(handle, path, columns, max_line_bytes)


def leading_comments(path):
    """Preserve native provenance comments, applying the same line bound."""
    with Path(path).open("rb") as handle:
        for _, line in bounded_lines(handle, path):
            if line.startswith("#"):
                yield line + "\n"
            elif line:
                return


def write_row(handle, fields, *, max_line_bytes=MAX_LINE_BYTES):
    """Write literal TSV without CSV quoting; bound merged output lines too."""
    fields = list(fields)
    size = 0
    for field in fields:
        raw = field.encode("utf-8")
        if b"\t" in raw or b"\n" in raw or CONTROL.search(raw):
            raise ValueError(f"Invalid delimiter/control character in featureCounts output: {handle.name}")
        size += len(raw) + 1  # One tab per field, with the last replaced by LF.
        if size > max_line_bytes:
            raise ValueError(f"featureCounts output line exceeds {max_line_bytes} bytes: {handle.name}")
    handle.write("\t".join(fields) + "\n")
