"""Paired fragment profiles honor the kit's reverse-stranded orientation."""
import csv
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gene_position_profile as profile


class FragmentProfileTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.gff = self.root / "features.gff"
        # Opposite-strand genes occupy the same interval. Directionality must
        # resolve them rather than report an ambiguous midpoint overlap.
        self.gff.write_text(
            "chr1\ttest\tCDS\t101\t300\t.\t+\t0\tlocus_tag=plus\n"
            "chr1\ttest\tCDS\t101\t300\t.\t-\t0\tlocus_tag=minus\n"
        )
        self.bam = self.root / "fragments.bam"
        self.prefix = self.root / "sample"

    def write_bam(self, pairs=True):
        with pysam.AlignmentFile(str(self.bam), "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}) as out:
            specs = [("plus", 163, 120, 220, 150), ("minus", 99, 120, 220, 150),
                     ("plus", 83, 220, 120, -150), ("minus", 147, 220, 120, -150)] if pairs else [("single", 16, 120, 0, 0)]
            for name, flag, start, mate_start, length in specs:
                read = pysam.AlignedSegment(out.header)
                read.query_name = name
                read.flag = flag
                read.reference_id = 0
                read.reference_start = start
                read.mapping_quality = 42
                read.cigarstring = "50M"
                read.query_sequence = "A" * 50
                if pairs:
                    read.next_reference_id = 0
                    read.next_reference_start = mate_start
                    read.template_length = length
                out.write(read)
        pysam.index(str(self.bam))

    def run_profile(self):
        args = ["profile", "--bam", str(self.bam), "--gff", str(self.gff),
                "--gene-id", "locus_tag", "--output-prefix", str(self.prefix),
                "--bins", "10", "--min-feature-reads", "1"]
        with patch.object(sys, "argv", args), patch.object(profile, "plot_profiles"):
            profile.main()

    def test_each_pair_counts_once_on_r2_sense_strand(self):
        self.write_bam()
        self.run_profile()
        with Path(str(self.prefix) + ".feature_read_positions.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 2)
        self.assertEqual({(row["fragment_id"], row["feature_id"]) for row in rows},
                         {("plus", "plus"), ("minus", "minus")})
        self.assertEqual({row["assignment_status"] for row in rows}, {"unique"})
        self.assertEqual({float(row["alignment_midpoint"]) for row in rows}, {195})

    def test_single_end_alignment_is_rejected(self):
        self.write_bam(pairs=False)
        with self.assertRaisesRegex(SystemExit, "unpaired alignment"):
            self.run_profile()


if __name__ == "__main__":
    unittest.main()
