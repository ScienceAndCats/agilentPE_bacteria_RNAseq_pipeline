import gzip
import json
from pathlib import Path
import tempfile
import unittest

import pysam

from mbc_processing import deduplicate_pairs, mbc_from_name, molecule_key, prepare_pairs


class MbcPreparationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        self.r1, self.r2 = self.base / "r1.fastq", self.base / "r2.fastq.gz"
        self.o1, self.o2 = self.base / "o1.fastq.gz", self.base / "o2.fastq"

    def inputs(self, names=("pair/1", "pair/2"), sequences=("ACGTAATTGCAA", "TGCGCTTAACC"), qualities=None):
        qualities = qualities or tuple("I" * len(s) for s in sequences)
        for path, name, sequence, quality in zip((self.r1, self.r2), names, sequences, qualities):
            opener = gzip.open if str(path).endswith(".gz") else open
            with opener(path, "wt") as handle:
                handle.write(f"@{name}\n{sequence}\n+\n{quality}\n")

    def test_default_strips_five_from_each_mate_without_deduplication(self):
        self.inputs()
        stats_path = self.base / "stats.json"
        result = prepare_pairs(self.r1, self.r2, self.o1, self.o2, stats_path=stats_path)
        self.assertEqual(gzip.open(self.o1, "rt").read(), "@pair\nATTGCAA\n+\nIIIIIII\n")
        self.assertEqual(self.o2.read_text(), "@pair\nTTAACC\n+\nIIIIII\n")
        self.assertEqual(result["pairs"], 1)
        self.assertEqual(json.loads(stats_path.read_text()), result)

    def test_count_extracts_three_bases_and_both_quality_strings(self):
        self.inputs(qualities=("!#$%&IIIIIII", "ABCDGIIIIII"))
        prepare_pairs(self.r1, self.r2, self.o1, self.o2, mode="count")
        name = self.o2.read_text().splitlines()[0][1:]
        self.assertEqual(mbc_from_name(name), ("ACG-TGC", "!#$ABC"))
        self.assertIn("\nATTGCAA\n", gzip.open(self.o1, "rt").read())
        self.assertEqual(gzip.open(self.o1, "rt").read().splitlines()[0], "@" + name)

    def test_already_trimmed_does_not_trim_again_and_count_fails(self):
        self.inputs()
        prepare_pairs(self.r1, self.r2, self.o1, self.o2, input_state="trimmed")
        self.assertIn("\nACGTAATTGCAA\n", gzip.open(self.o1, "rt").read())
        with self.assertRaisesRegex(ValueError, "requires raw"):
            prepare_pairs(self.r1, self.r2, self.o1, self.o2, mode="count", input_state="trimmed")

    def test_ambiguous_mbc_rejected_only_for_count(self):
        self.inputs(sequences=("NAGCCATTGCAA", "TGCGCTTAACC"))
        with self.assertRaisesRegex(ValueError, "ambiguous MBC"):
            prepare_pairs(self.r1, self.r2, self.o1, self.o2, mode="count")
        prepare_pairs(self.r1, self.r2, self.o1, self.o2)

    def test_fastq_pair_validation(self):
        for names in (("first/1", "second/2"), ("pair/2", "pair/1"), ("pair 2:N:0:A", "pair 2:N:0:A")):
            with self.subTest(names=names):
                self.inputs(names=names)
                with self.assertRaises(ValueError):
                    prepare_pairs(self.r1, self.r2, self.o1, self.o2)
        self.inputs()
        with gzip.open(self.r2, "wt") as handle:
            handle.write("")
        with self.assertRaisesRegex(ValueError, "different numbers"):
            prepare_pairs(self.r1, self.r2, self.o1, self.o2)

    def test_malformed_fastq_and_too_short_are_rejected(self):
        self.inputs(qualities=("I", "I" * 11))
        with self.assertRaisesRegex(ValueError, "Malformed FASTQ"):
            prepare_pairs(self.r1, self.r2, self.o1, self.o2)
        self.inputs(sequences=("ACGTA", "TGCGCTTAACC"))
        with self.assertRaisesRegex(ValueError, "no insert"):
            prepare_pairs(self.r1, self.r2, self.o1, self.o2)

    def test_overwriting_inputs_and_double_extraction_are_rejected(self):
        self.inputs()
        with self.assertRaisesRegex(ValueError, "distinct"):
            prepare_pairs(self.r1, self.r2, self.r1, self.o2)
        self.inputs(names=("pair:MBC:ACG-TGC:494949494949", "pair:MBC:ACG-TGC:494949494949"))
        with self.assertRaisesRegex(ValueError, "already contains"):
            prepare_pairs(self.r1, self.r2, self.o1, self.o2)


HEADER = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 10000}, {"SN": "chr2", "LN": 10000}]}


def pair(name, barcode="ACG-TGC", start1=100, start2=200, reference=0, reverse=False, unmapped=False):
    qname = f"{name}:MBC:{barcode}:4949494a4a4a"
    records = []
    for first in (True, False):
        record = pysam.AlignedSegment()
        record.query_name = qname
        record.query_sequence = "A" * 30
        record.query_qualities = pysam.qualitystring_to_array("I" * 30)
        record.flag = (83 if first else 163) if reverse else (99 if first else 147)
        record.reference_id = reference
        record.reference_start = start1 if first else start2
        record.next_reference_id = reference
        record.next_reference_start = start2 if first else start1
        record.cigarstring = "30M"
        record.mapping_quality = 42
        record.template_length = 130 if first else -130
        if unmapped:
            record.flag = 77 if first else 141
            record.reference_id = record.next_reference_id = -1
            record.reference_start = record.next_reference_start = -1
            record.cigarstring = None
        records.append(record)
    return records


class MbcDeduplicationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)

    def run_dedup(self, records, stem="sample"):
        source, output = self.base / f"{stem}.sam", self.base / f"{stem}.bam"
        with pysam.AlignmentFile(str(source), "w", header=HEADER) as handle:
            for record in records:
                handle.write(record)
        result = deduplicate_pairs(source, output, self.base / f"{stem}.json")
        with pysam.AlignmentFile(str(output), "rb") as handle:
            kept = list(handle.fetch(until_eof=True))
        return result, kept

    def test_exact_dual_barcode_and_both_endpoints_define_molecules(self):
        groups = [pair("duplicate1"), pair("duplicate2"), pair("other_barcode", "ACG-TGA"),
                  pair("other_start", start1=101), pair("other_end", start2=201),
                  pair("other_reference", reference=1), pair("other_strand", reverse=True),
                  pair("unmapped", unmapped=True)]
        # Coordinate-like interleaving: pair mates do not have to be adjacent.
        result, kept = self.run_dedup([r for mate in (0, 1) for group in reversed(groups) for r in [group[mate]]])
        self.assertEqual(result["input_pairs"], 8)
        self.assertEqual(result["mapped_pairs"], 7)
        self.assertEqual(result["unmapped_pairs"], 1)
        self.assertEqual(result["molecules"], 6)
        self.assertEqual(result["duplicate_pairs"], 1)
        self.assertEqual(len(kept), 12)
        self.assertTrue(all(r.has_tag("RX") and r.has_tag("QX") and r.has_tag("MI") for r in kept))
        self.assertEqual(kept[0].get_tag("RX"), "ACG-TGC")
        self.assertEqual(kept[0].get_tag("QX"), "III JJJ")

    def test_independent_samples_do_not_collapse_each_other(self):
        for stem in ("sample1", "sample2"):
            result, kept = self.run_dedup(pair("same"), stem)
            self.assertEqual(result["molecules"], 1)
            self.assertEqual(len(kept), 2)

    def test_missing_mate_and_missing_mbc_fail(self):
        with self.assertRaisesRegex(ValueError, "exactly one primary"):
            self.run_dedup(pair("orphan")[:1], "orphan")
        records = pair("untagged")
        for record in records:
            record.query_name = "untagged"
        with self.assertRaisesRegex(ValueError, "lacks valid dual MBC"):
            self.run_dedup(records, "untagged")

    def test_clipping_uses_unclipped_five_prime_coordinates(self):
        original, equivalent = pair("original"), pair("clipped", start1=105, start2=200)
        equivalent[0].cigarstring = "5S25M"
        # Read 2 is reverse, so 5′ soft clipping belongs at the right edge.
        equivalent[1].cigarstring = "25M5S"
        self.assertEqual(molecule_key(*original, "ACG-TGC"), molecule_key(*equivalent, "ACG-TGC"))
        result, kept = self.run_dedup(original + equivalent)
        self.assertEqual(result["molecules"], 1)
        self.assertEqual(len(kept), 2)

    def test_secondary_records_do_not_create_extra_molecules(self):
        records = pair("primary")
        secondary = pair("primary")[0]
        secondary.flag |= 256
        result, kept = self.run_dedup(records + [secondary])
        self.assertEqual(result["excluded_nonprimary_alignments"], 1)
        self.assertEqual(result["molecules"], 1)


if __name__ == "__main__":
    unittest.main()
