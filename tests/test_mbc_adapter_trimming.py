"""Exercise the actual pipeline preprocessing block with the installed Cutadapt."""
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest

from mbc_processing import mbc_from_name


REPOSITORY = Path(__file__).resolve().parents[1]
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
INSERT = "GTCACGATTCGATCGTACGTTGCACTAGCGATTCGACATG"
LEADER_R1, LEADER_R2 = "ACGTA", "TGCGC"


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def fastq_records(path):
    lines = path.read_text().splitlines()
    return [(lines[index][1:], lines[index + 1], lines[index + 3]) for index in range(0, len(lines), 4)]


@unittest.skipUnless(shutil.which("cutadapt"), "real Cutadapt is required")
class MbcAdapterTrimmingTests(unittest.TestCase):
    def run_preprocessing(self, sequences, *, mode="strip", input_state="raw"):
        temporary = tempfile.TemporaryDirectory(prefix="mbc adapter ")
        self.addCleanup(temporary.cleanup)
        base = Path(temporary.name)
        input_paths = [base / "reads R1.fastq", base / "reads R2.fastq"]
        for mate, path in enumerate(input_paths):
            with path.open("w") as handle:
                for number, pair in enumerate(sequences):
                    sequence = pair[mate]
                    handle.write(f"@pair{number}/{mate + 1}\n{sequence}\n+\n{'I' * len(sequence)}\n")
        script = (REPOSITORY / "map_bacteria_with_decoys.sh").read_text()
        match = re.search(r'^mkdir -p "\$OUTPUT_DIR/\.mbc_fastq".*?^done$', script, re.MULTILINE | re.DOTALL)
        self.assertIsNotNone(match, "Could not locate the pipeline MBC/Cutadapt processing block")
        # Execute the real shell block, including the actual adapter patterns,
        # rather than recreating a similar trimming command inside the test.
        assignments = {
            "OUTPUT_DIR": str(base), "SCRIPT_DIR": str(REPOSITORY), "MBC_MODE": mode,
            "MBC_INPUT": input_state, "THREADS": "1", "MIN_READ_LENGTH": "22",
        }
        commands = ["set -euo pipefail", f"source {shlex.quote(str(REPOSITORY / 'config.env'))}"]
        commands.extend(f"{key}={shlex.quote(value)}" for key, value in assignments.items())
        commands += ["sample_names=(sample)",
                     f"read1_filenames=({shlex.quote(str(input_paths[0]))})",
                     f"read2_filenames=({shlex.quote(str(input_paths[1]))})", match.group()]
        environment = os.environ.copy()
        environment["PATH"] = str(Path(sys.executable).parent) + os.pathsep + environment.get("PATH", "")
        result = subprocess.run(["bash", "-c", "\n".join(commands)], text=True, capture_output=True, env=environment)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        outputs = [fastq_records(base / f"sample_R{mate}_22bp.trim.fastq") for mate in (1, 2)]
        stats = json.loads((base / "sample.mbc_prepare.json").read_text())
        return outputs, stats

    def raw_with_readthrough(self, *, adapter_bases=None):
        adapter1 = ADAPTER_R1 if adapter_bases is None else ADAPTER_R1[:adapter_bases]
        adapter2 = ADAPTER_R2 if adapter_bases is None else ADAPTER_R2[:adapter_bases]
        return (LEADER_R1 + INSERT + reverse_complement(LEADER_R2) + adapter1,
                LEADER_R2 + reverse_complement(INSERT) + reverse_complement(LEADER_R1) + adapter2)

    def test_strip_removes_both_leaders_and_opposite_mbc_dark_readthrough(self):
        outputs, stats = self.run_preprocessing([self.raw_with_readthrough(), self.raw_with_readthrough()])
        for mate, expected in zip(outputs, (INSERT, reverse_complement(INSERT))):
            self.assertEqual(len(mate), 2, "Strip mode must keep duplicate read pairs")
            self.assertEqual([record[1] for record in mate], [expected, expected])
            self.assertTrue(all(len(record[2]) == len(expected) for record in mate))
            self.assertTrue(all(":MBC:" not in record[0] for record in mate))
        self.assertEqual(stats["pairs"], 2)

    def test_count_preserves_dual_barcodes_while_trimming_readthrough(self):
        outputs, stats = self.run_preprocessing([self.raw_with_readthrough(), self.raw_with_readthrough()], mode="count")
        self.assertEqual([len(mate) for mate in outputs], [2, 2])
        for first, second in zip(*outputs):
            self.assertEqual(first[0], second[0])
            self.assertEqual(mbc_from_name(first[0]), ("ACG-TGC", "IIIIII"))
            self.assertEqual((first[1], second[1]), (INSERT, reverse_complement(INSERT)))
        self.assertEqual(stats["mode"], "count")
        self.assertEqual(stats["pairs"], 2)

    def test_eight_adapter_bases_suffice_to_remove_five_opposite_bases(self):
        outputs, _ = self.run_preprocessing([self.raw_with_readthrough(adapter_bases=8)])
        self.assertEqual((outputs[0][0][1], outputs[1][0][1]), (INSERT, reverse_complement(INSERT)))

    def test_already_trimmed_input_preserves_insert_leading_bases(self):
        pair = (INSERT + reverse_complement(LEADER_R2) + ADAPTER_R1,
                reverse_complement(INSERT) + reverse_complement(LEADER_R1) + ADAPTER_R2)
        outputs, stats = self.run_preprocessing([pair], input_state="trimmed")
        self.assertEqual((outputs[0][0][1], outputs[1][0][1]), (INSERT, reverse_complement(INSERT)))
        self.assertEqual(stats["leading_bases_removed_per_mate"], 0)

    def test_opposite_only_tail_without_adapter_evidence_is_conservatively_retained(self):
        # A five-base opposite-end prefix alone cannot distinguish read-through
        # from genuine insert sequence. The current matcher requires adapter
        # evidence rather than deleting potentially biological five-base matches.
        for length in (2, 5):
            with self.subTest(tail_bases=length):
                pair = (LEADER_R1 + INSERT + reverse_complement(LEADER_R2)[:length],
                        LEADER_R2 + reverse_complement(INSERT) + reverse_complement(LEADER_R1)[:length])
                outputs, _ = self.run_preprocessing([pair])
                self.assertEqual(outputs[0][0][1], pair[0][5:])
                self.assertEqual(outputs[1][0][1], pair[1][5:])


if __name__ == "__main__":
    unittest.main()
