import csv
import io
from pathlib import Path
import tempfile
import unittest

from featurecounts_io import (
    ANNOTATION_COLUMNS, MAX_COUNT, MAX_LINE_BYTES, MAX_SUMMARY_ROWS,
    bounded_lines, open_table, write_row,
)


HEADER = "\t".join((*ANNOTATION_COLUMNS, "sample.sam")) + "\n"
ROW = "gene\tchr1\t1\t20\t+\t20\t3\n"


class FeatureCountsReaderTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.path = Path(temporary.name) / "counts.txt"

    def rows(self, data, *, summary=False, **kwargs):
        self.path.write_bytes(data if isinstance(data, bytes) else data.encode("utf-8"))
        with open_table(self.path, ("Status",) if summary else ANNOTATION_COLUMNS, **kwargs) as rows:
            return list(rows)

    def test_native_long_annotation_and_literal_quotes_without_changing_csv_limit(self):
        old_limit = csv.field_size_limit()
        annotation = ";".join(["chr1"] * 40000)
        data = HEADER + f'gene"with,quotes\t{annotation}\t1\t20\t+\t20\t3\n'
        rows = self.rows(data)
        self.assertEqual(rows[1][0], 'gene"with,quotes')
        self.assertEqual(rows[1][1], annotation)
        self.assertEqual(rows[1][-1], "3")
        self.assertEqual(csv.field_size_limit(), old_limit)

    def test_bounded_read_fails_without_reading_the_rest_of_an_overlong_line(self):
        class TrackedFile(io.BytesIO):
            requests = []
            def readline(self, size=-1):
                self.requests.append(size)
                return super().readline(size)
        handle = TrackedFile(b"x" * 1000)
        with self.assertRaisesRegex(ValueError, "exceeds 64 bytes at test:1"):
            list(bounded_lines(handle, "test", max_line_bytes=64))
        self.assertEqual(handle.requests, [65])
        self.assertEqual(handle.tell(), 65)

    def test_default_cap_covers_comments_headers_and_data_even_without_newline(self):
        for prefix in (b"#", b"", HEADER.encode() + b"gene\t"):
            with self.subTest(prefix=prefix):
                with self.assertRaisesRegex(ValueError, f"exceeds {MAX_LINE_BYTES} bytes"):
                    self.rows(prefix + b"A" * (MAX_LINE_BYTES + 1))

    def test_exact_byte_limit_crlf_and_utf8(self):
        line = "é\tvalue\r\n".encode("utf-8")
        self.assertEqual(list(bounded_lines(io.BytesIO(line), "test", max_line_bytes=len(line))),
                         [(1, "é\tvalue")])
        with self.assertRaisesRegex(ValueError, "exceeds"):
            list(bounded_lines(io.BytesIO(line), "test", max_line_bytes=len(line) - 1))
        self.assertEqual(self.rows((HEADER + ROW).replace("\n", "\r\n"))[1][-1], "3")

    def test_bad_encoding_and_control_characters_have_path_and_line(self):
        for bad in (b"\xff", b"\x00", b"\r", b"\x1b"):
            with self.subTest(bad=bad):
                with self.assertRaisesRegex(ValueError, rf"{self.path}:2"):
                    self.rows(HEADER.encode() + b"gene" + bad + ROW[4:].encode())

    def test_header_and_row_shape_are_validated(self):
        cases = ["", "# comment only\n", HEADER.replace("Geneid", "Other"),
                 HEADER.replace("sample.sam", ""),
                 HEADER.rstrip("\n") + "\tsample.sam\n",
                 HEADER + ROW.rstrip("\n") + "\textra\n",
                 HEADER + "\t".join(ROW.strip().split("\t")[:-1]) + "\n",
                 HEADER + ROW.replace("chr1", "")]
        for data in cases:
            with self.subTest(data=data[:100]), self.assertRaises(ValueError):
                self.rows(data)

    def test_only_nonnegative_unsigned_64_bit_counts_and_lengths_are_allowed(self):
        self.assertEqual(self.rows(HEADER + ROW.replace("\t3\n", f"\t{MAX_COUNT}\n"))[1][-1], str(MAX_COUNT))
        for value in ("-1", "1.5", "NaN", "1e9", "", "٣", str(MAX_COUNT + 1), "9" * 10000):
            with self.subTest(value=value[:30]), self.assertRaisesRegex(ValueError, "Invalid featureCounts count"):
                self.rows(HEADER + ROW.replace("\t3\n", f"\t{value}\n"))
        with self.assertRaisesRegex(ValueError, "feature length"):
            self.rows(HEADER + ROW.replace("\t20\t3", "\tbad\t3"))

    def test_summary_requires_assigned_and_unique_bounded_statuses(self):
        header = "Status\tsample.sam\n"
        self.assertEqual(self.rows(header + "Assigned\t3\nUnassigned_Unmapped\t0\n", summary=True)[1], ["Assigned", "3"])
        for body in ("Unassigned_Unmapped\t0\n", "Assigned\t3\nAssigned\t0\n",
                     "Assigned\t-1\n", "Assigned\t3\textra\n", "bad status\t0\n"):
            with self.subTest(body=body), self.assertRaises(ValueError):
                self.rows(header + body, summary=True)
        body = "Assigned\t0\n" + "".join(f"Unassigned_{i}\t0\n" for i in range(MAX_SUMMARY_ROWS))
        with self.assertRaisesRegex(ValueError, "Too many"):
            self.rows(header + body, summary=True)

    def test_writer_preserves_quotes_and_checks_output_size_before_writing(self):
        with self.path.open("w", encoding="utf-8") as handle:
            write_row(handle, ['gene"quoted,', "3"])
        self.assertEqual(self.path.read_text(), 'gene"quoted,\t3\n')
        with self.path.open("w", encoding="utf-8") as handle:
            with self.assertRaisesRegex(ValueError, "output line exceeds"):
                write_row(handle, ["é" * 10, "1"], max_line_bytes=20)
            for invalid in ("bad\tvalue", "bad\nvalue", "bad\rvalue", "bad\x00value"):
                with self.assertRaises(ValueError):
                    write_row(handle, [invalid, "1"])
        self.assertEqual(self.path.read_bytes(), b"")


if __name__ == "__main__":
    unittest.main()
