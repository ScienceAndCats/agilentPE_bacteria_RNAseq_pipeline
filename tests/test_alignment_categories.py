import json
from pathlib import Path
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from alignment_categories import AnnotationIndex, classify_alignments


class CategoryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.gff = self.root / "features.gff"
        self.gff.write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t1\t1000\t.\t+\t.\tID=gene;gene=gene\n"
            "chr1\ttest\tmRNA\t1\t1000\t.\t+\t.\tID=transcript;Parent=gene\n"
            "chr1\ttest\texon\t1\t200\t.\t+\t.\tParent=transcript\n"
            "chr1\ttest\texon\t801\t1000\t.\t+\t.\tParent=transcript\n"
            "chr1\ttest\tCDS\t101\t200\t.\t+\t0\tlocus_tag=cds1\n"
            "chr1\ttest\tCDS\t801\t900\t.\t+\t0\tlocus_tag=cds1\n"
            "chr1\ttest\trRNA\t150\t200\t.\t+\t.\tlocus_tag=rrna1\n"
            "chr1\ttest\ttRNA\t1201\t1300\t.\t-\t.\tlocus_tag=trna1\n"
        )
        self.sam = self.root / "reads.sam"

    def write_reads(self, specs):
        with pysam.AlignmentFile(str(self.sam), "w", header={"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 2000}]}) as handle:
            for name, start, cigar, flag in specs:
                read = pysam.AlignedSegment(handle.header)
                read.query_name = name
                read.flag = flag
                if not flag & 4:
                    read.reference_id = 0
                    read.reference_start = start
                    read.mapping_quality = 42
                    read.cigarstring = cigar
                    if flag & 1:
                        read.next_reference_id = 0
                        read.next_reference_start = start
                        read.template_length = 50 if flag & 64 else -50
                    read.query_sequence = "A" * read.infer_query_length()
                else:
                    read.query_sequence = "A" * 50
                handle.write(read)

    def write_pairs(self, specs):
        """Place both mates at the selected locus on the library's sense strand."""
        reads = []
        for name, start, cigar, strand in specs:
            flags = (83, 163) if strand == "+" else (99, 147)
            reads.extend((name, start, cigar, flag) for flag in flags)
        self.write_reads(reads)

    def test_categories_overlap_priority_and_gene_intron(self):
        self.write_pairs([
            ("cds", 100, "20M", "+"),
            ("rrna_priority", 150, "20M", "+"),
            ("utr", 20, "20M", "+"),
            ("intron", 400, "20M", "+"),
            ("unannotated", 1500, "20M", "+"),
            ("trna", 1200, "20M", "-"),
        ])
        result = classify_alignments(self.sam, self.gff)
        self.assertEqual(result["aligned"], 6)
        self.assertEqual(result["categories"], {"rRNA": 1, "CDS": 1, "other_RNA": 2, "intergenic_unannotated": 2})
        self.assertEqual(result["cds_features_with_reads"], ["cds1"])

    def test_spliced_alignment_does_not_count_unsequenced_gap(self):
        self.write_pairs([("spliced", 10, "20M1170N20M", "+")])
        # Neither block touches CDS/rRNA; the intervening skipped span touches both.
        result = classify_alignments(self.sam, self.gff)
        self.assertEqual(result["categories"]["other_RNA"], 1)
        self.assertEqual(result["categories"]["rRNA"], 0)
        self.assertEqual(result["categories"]["CDS"], 0)

    def test_pair_counted_once_and_fragment_gap_excluded(self):
        self.write_reads([
            ("pair", 10, "20M", 83),
            ("pair", 950, "20M", 163),
            ("discordant", 100, "20M", 65),
            ("discordant", 150, "20M", 129),
        ])
        result = classify_alignments(self.sam, self.gff)
        self.assertEqual(result["aligned"], 1)
        self.assertEqual(result["categories"]["other_RNA"], 1)
        self.assertEqual(result["cds_features_with_reads"], [])

    def test_r1_antisense_and_r2_sense_select_transcript_strand(self):
        self.write_pairs([("sense_transcript", 100, "20M", "+"),
                          ("opposite_transcript", 100, "20M", "-")])
        result = classify_alignments(self.sam, self.gff)
        self.assertEqual(result["categories"]["CDS"], 1)
        self.assertEqual(result["categories"]["intergenic_unannotated"], 1)

    def test_unpaired_primary_alignment_is_rejected(self):
        self.write_reads([("single", 100, "20M", 16)])
        with self.assertRaisesRegex(ValueError, "Unpaired alignment"):
            classify_alignments(self.sam, self.gff)

    def test_nonprimary_and_unmapped_alignments_are_ignored(self):
        self.write_reads([("unmapped", 0, "", 4),
                          ("secondary", 100, "20M", 256),
                          ("supplementary", 100, "20M", 2048)])
        self.assertEqual(classify_alignments(self.sam, self.gff)["aligned"], 0)

    def test_rna_transcript_utr_without_exon_records(self):
        self.gff.write_text("chr1\ttest\tmRNA\t1\t1000\t.\t+\t.\tID=transcript\nchr1\ttest\tCDS\t101\t200\t.\t+\t0\tlocus_tag=cds1\n")
        self.write_pairs([("utr", 20, "20M", "+")])
        self.assertEqual(classify_alignments(self.sam, self.gff, reference_kind="rna")["categories"]["other_RNA"], 1)
        self.assertEqual(classify_alignments(self.sam, self.gff, reference_kind="dna")["categories"]["intergenic_unannotated"], 1)

    def test_missing_adjacent_mate_fails(self):
        self.write_reads([("pair", 10, "20M", 99)])
        with self.assertRaisesRegex(ValueError, "adjacent primary mates"):
            classify_alignments(self.sam, self.gff)

    def test_explicit_utrs_count_as_other_rna(self):
        self.gff.write_text("chr1\ttest\tfive_prime_UTR\t1\t99\t.\t+\t.\tID=utr1\nchr1\ttest\tthree_prime_UTR\t901\t1000\t.\t+\t.\tID=utr2\n")
        self.write_pairs([("fiveprime", 10, "20M", "+"), ("threeprime", 950, "20M", "+")])
        self.assertEqual(classify_alignments(self.sam, self.gff)["categories"]["other_RNA"], 2)

    def test_mismatched_reference_names_fail_even_without_reads(self):
        self.gff.write_text("different_genome\ttest\tCDS\t1\t100\t.\t+\t.\tID=cds\n")
        self.write_reads([])
        with self.assertRaisesRegex(ValueError, "reference names.*absent"):
            classify_alignments(self.sam, self.gff)

    def test_missing_relevant_annotation_fails(self):
        self.gff.write_text("chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=gene\n")
        self.write_reads([])
        with self.assertRaisesRegex(ValueError, "No relevant"):
            classify_alignments(self.sam, self.gff)

    def test_annotation_beyond_reference_end_fails(self):
        self.gff.write_text("chr1\ttest\tCDS\t1901\t2100\t.\t+\t.\tID=cds\n")
        self.write_reads([])
        with self.assertRaisesRegex(ValueError, "coordinates exceed"):
            classify_alignments(self.sam, self.gff)

    def test_interval_boundaries_and_unknown_strand(self):
        self.gff.write_text("chr1\ttest\tCDS\t100\t200\t.\t.\t0\tlocus_tag=cds\n")
        index = AnnotationIndex(self.gff)
        self.assertEqual(index.classify("chr1", [(0, 99)], "+")[0], "intergenic_unannotated")
        self.assertEqual(index.classify("chr1", [(200, 250)], "-")[0], "intergenic_unannotated")
        self.assertEqual(index.classify("chr1", [(99, 100)], "-")[0], "CDS")


if __name__ == "__main__":
    unittest.main()
