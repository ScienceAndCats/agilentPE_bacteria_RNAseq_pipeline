"""SureSelect Max inline MBC extraction and exact paired-fragment deduplication.

The kit has three barcode bases and two dark bases at the start of each mate
(Agilent G9689-90000, p. 31). Deduplication estimates captured fragments; it does
not perform barcode error correction, build consensus sequences, or recover
barcode collisions. All operations are sample-local.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
import gzip
import hashlib
from itertools import groupby, zip_longest
import json
from pathlib import Path
import re
import sqlite3
import tempfile

import pysam


MBC_MARKER = ":MBC:"
MBC_NAME = re.compile(r"^(?P<name>.+):MBC:(?P<barcode>[ACGT]{3}-[ACGT]{3}):(?P<quality>[0-9a-f]{12})$")


def _open_fastq(path, mode):
    return gzip.open(path, mode + "t") if str(path).endswith(".gz") else open(path, mode)


def _records(handle, path):
    number = 0
    while True:
        header = handle.readline()
        if not header:
            return
        number += 1
        sequence, plus, quality = (handle.readline().rstrip("\r\n") for _ in range(3))
        header = header.rstrip("\r\n")
        if (not header.startswith("@") or not header[1:].split() or not plus.startswith("+")
                or not sequence or len(sequence) != len(quality)
                or any(not 33 <= ord(c) <= 126 for c in quality)):
            raise ValueError(f"Malformed FASTQ record {number} in {path}")
        yield header[1:], sequence, quality


def _identity(header, mate):
    fields = header.split(maxsplit=1)
    name = fields[0]
    if name.endswith(("/1", "/2")):
        if not name.endswith(f"/{mate}"):
            raise ValueError(f"Read {mate} has incorrect mate suffix: {name}")
        name = name[:-2]
    comment = fields[1] if len(fields) > 1 else ""
    if re.match(r"^[12]:", comment) and not comment.startswith(f"{mate}:"):
        raise ValueError(f"Read {mate} has incorrect mate field: {header}")
    return name, comment


def _write_stats(path, stats):
    if path is not None:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile("w", dir=destination.parent, delete=False) as handle:
            json.dump(stats, handle, indent=2)
            handle.write("\n")
            temporary = Path(handle.name)
        temporary.replace(destination)


def prepare_pairs(read1, read2, out1, out2, mode="strip", input_state="raw", stats_path=None):
    """Validate mates; trim 5′ kit bases, optionally preserving both MBCs in QNAME.

    ``trimmed`` means upstream demultiplexing already removed all five leading
    bases; it is accepted only for strip mode to prevent irrecoverable UMI loss.
    """
    if mode not in {"strip", "count"} or input_state not in {"raw", "trimmed"}:
        raise ValueError("MBC mode must be strip/count and input must be raw/trimmed")
    if mode == "count" and input_state != "raw":
        raise ValueError("MBC count mode requires raw reads with the inline barcodes intact")
    inputs = {Path(read1).resolve(), Path(read2).resolve()}
    outputs = {Path(out1).resolve(), Path(out2).resolve()}
    if len(inputs) != 2 or len(outputs) != 2 or inputs & outputs:
        raise ValueError("MBC inputs and outputs must be distinct files")
    stats = {"mode": mode, "input": input_state, "pairs": 0,
             "leading_bases_removed_per_mate": 5 if input_state == "raw" else 0,
             "barcode_bases_per_mate": 3, "dark_bases_per_mate": 2}
    with ExitStack() as stack:
        readers = [stack.enter_context(_open_fastq(p, "r")) for p in (read1, read2)]
        writers = [stack.enter_context(_open_fastq(p, "w")) for p in (out1, out2)]
        for pair in zip_longest(*[_records(h, p) for h, p in zip(readers, (read1, read2))]):
            if any(record is None for record in pair):
                raise ValueError("Paired FASTQs contain different numbers of records")
            identities = [_identity(record[0], mate) for mate, record in enumerate(pair, 1)]
            if identities[0][0] != identities[1][0]:
                raise ValueError(f"Mismatched paired FASTQ identifiers: {pair[0][0]!r}, {pair[1][0]!r}")
            name = identities[0][0]
            if MBC_MARKER in name:
                raise ValueError(f"Read name already contains MBC metadata: {name}")
            if input_state == "raw" and any(len(r[1]) <= 5 for r in pair):
                raise ValueError(f"Read pair {name} has no insert after removing five inline bases")
            if mode == "count":
                barcodes = [r[1][:3].upper() for r in pair]
                if any(not re.fullmatch(r"[ACGT]{3}", barcode) for barcode in barcodes):
                    raise ValueError(f"Read pair {name} has an ambiguous MBC; count mode requires A/C/G/T barcodes")
                quality = "".join(r[2][:3] for r in pair).encode("ascii").hex()
                name += f"{MBC_MARKER}{barcodes[0]}-{barcodes[1]}:{quality}"
            if len(name.encode("ascii")) > 254:
                raise ValueError(f"Read name exceeds the SAM QNAME length limit after MBC extraction: {name}")
            for writer, record, identity in zip(writers, pair, identities):
                sequence, quality = record[1:]
                if input_state == "raw":
                    sequence, quality = sequence[5:], quality[5:]
                comment = f" {identity[1]}" if identity[1] else ""
                writer.write(f"@{name}{comment}\n{sequence}\n+\n{quality}\n")
            stats["pairs"] += 1
    _write_stats(stats_path, stats)
    return stats


def mbc_from_name(name):
    match = MBC_NAME.fullmatch(name)
    if not match:
        raise ValueError(f"Aligned read lacks valid dual MBC metadata: {name}")
    qualities = bytes.fromhex(match["quality"]).decode("ascii")
    if any(not 33 <= ord(c) <= 126 for c in qualities):
        raise ValueError(f"Invalid MBC base qualities in aligned read: {name}")
    return match["barcode"], qualities


def _five_prime(record):
    """Zero-based unclipped 5′ coordinate, allowing negative clipped positions."""
    cigar = record.cigartuples or []
    if record.is_reverse:
        clipped = 0
        for operation, length in reversed(cigar):
            if operation not in (4, 5):
                break
            clipped += length
        return record.reference_end - 1 + clipped
    clipped = 0
    for operation, length in cigar:
        if operation not in (4, 5):
            break
        clipped += length
    return record.reference_start - clipped


def molecule_key(first, second, barcode):
    return (barcode, first.reference_id, _five_prime(first), first.is_reverse,
            second.reference_id, _five_prime(second), second.is_reverse)


def deduplicate_pairs(input_path, output_path, stats_path=None):
    """Write one mapped primary pair per exact MBC/fragment identity to BAM.

    Input need not be sorted. Name sorting and molecule keys use temporary disk
    storage, keeping memory bounded. The first pair in query-name order represents
    each family; sequences are not combined into a consensus. Unmapped pairs and
    nonprimary alignments are excluded from the molecule-only BAM.
    """
    if Path(input_path).resolve() == Path(output_path).resolve():
        raise ValueError("Deduplication input and output must be distinct files")
    stats = {"method": "exact_dual_mbc_and_fragment_endpoints", "input_pairs": 0,
             "mapped_pairs": 0, "unmapped_pairs": 0, "molecules": 0, "duplicate_pairs": 0,
             "excluded_nonprimary_alignments": 0}
    with tempfile.TemporaryDirectory(prefix="mbc-dedup-", dir=Path(output_path).resolve().parent) as directory:
        sorted_path = str(Path(directory) / "name_sorted.bam")
        pysam.sort("-n", "-m", "128M", "-T", str(Path(directory) / "sort"), "-o", sorted_path, str(input_path))
        database = sqlite3.connect(Path(directory) / "molecules.sqlite")
        database.execute("CREATE TABLE seen (identity TEXT PRIMARY KEY) WITHOUT ROWID")
        try:
            with pysam.AlignmentFile(sorted_path, "rb") as source:
                header = source.header.to_dict()
                header.setdefault("HD", {})["SO"] = "queryname"
                with pysam.AlignmentFile(str(output_path), "wb", header=header) as destination:
                    for name, alignments in groupby(source.fetch(until_eof=True), key=lambda r: r.query_name):
                        primary = []
                        for record in alignments:
                            if record.is_secondary or record.is_supplementary:
                                stats["excluded_nonprimary_alignments"] += 1
                            else:
                                primary.append(record)
                        if not primary:
                            continue
                        if (len(primary) != 2 or any(not r.is_paired for r in primary)
                                or sum(r.is_read1 for r in primary) != 1
                                or sum(r.is_read2 for r in primary) != 1):
                            raise ValueError(f"Expected exactly one primary R1/R2 pair for {name}")
                        first, second = sorted(primary, key=lambda r: not r.is_read1)
                        barcode, qualities = mbc_from_name(name)
                        stats["input_pairs"] += 1
                        if any(r.is_unmapped for r in primary):
                            stats["unmapped_pairs"] += 1
                            continue
                        if any(not r.is_proper_pair for r in primary):
                            raise ValueError(f"Molecule counting requires concordantly mapped pairs: {name}")
                        stats["mapped_pairs"] += 1
                        key = json.dumps(molecule_key(first, second, barcode), separators=(",", ":"))
                        inserted = database.execute("INSERT OR IGNORE INTO seen VALUES (?)", (key,)).rowcount
                        if not inserted:
                            stats["duplicate_pairs"] += 1
                            continue
                        stats["molecules"] += 1
                        identifier = hashlib.sha256(key.encode()).hexdigest()[:24]
                        for record in (first, second):
                            record.set_tag("RX", barcode, value_type="Z")
                            record.set_tag("QX", qualities[:3] + " " + qualities[3:], value_type="Z")
                            record.set_tag("MI", identifier, value_type="Z")
                            destination.write(record)
        finally:
            database.close()
    _write_stats(stats_path, stats)
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare", help="Remove kit prefixes, preserving MBCs in count mode")
    for name in ("read1", "read2", "out1", "out2"):
        prepare.add_argument("--" + name, required=True)
    prepare.add_argument("--mode", choices=("strip", "count"), default="strip")
    prepare.add_argument("--input", choices=("raw", "trimmed"), default="raw", dest="input_state")
    prepare.add_argument("--stats", dest="stats_path")
    dedup = commands.add_parser("deduplicate", help="Write exact-MBC deduplicated mapped paired BAM")
    dedup.add_argument("--input", required=True, dest="input_path")
    dedup.add_argument("--output", required=True, dest="output_path")
    dedup.add_argument("--stats", dest="stats_path")
    args = vars(parser.parse_args())
    command = args.pop("command")
    try:
        result = prepare_pairs(**args) if command == "prepare" else deduplicate_pairs(**args)
    except (ValueError, OSError, pysam.SamtoolsError) as error:
        parser.exit(1, f"ERROR: {error}\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
