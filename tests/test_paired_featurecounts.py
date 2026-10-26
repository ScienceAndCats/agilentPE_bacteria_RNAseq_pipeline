"""Real featureCounts checks for paired libraries with empty residual stages."""
import csv
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import pysam

from paired_featurecounts import add_zero_columns, count_features
from featurecounts_io import ANNOTATION_COLUMNS, open_table


class FeatureCountsPostprocessingTests(unittest.TestCase):
    def test_long_fields_and_quotes_survive_zero_column_insertion(self):
        limit = csv.field_size_limit()
        with tempfile.TemporaryDirectory() as temporary:
            source, output = Path(temporary) / 'source.txt', Path(temporary) / 'final.txt'
            annotation = ';'.join(['chr1'] * 40000)
            source.write_text('# provenance\n' + '\t'.join((*ANNOTATION_COLUMNS, 'reads.sam')) +
                              f'\ngene"quoted,\t{annotation}\t1\t20\t+\t20\t7\n')
            add_zero_columns(source, output, ['empty.bam', 'reads.sam'], ['reads.sam'], ANNOTATION_COLUMNS)
            with open_table(output) as rows:
                header, row = list(rows)
            self.assertEqual(header[6:], ['empty.bam', 'reads.sam'])
            self.assertEqual(row, ['gene"quoted,', annotation, '1', '20', '+', '20', '0', '7'])
            self.assertTrue(output.read_text().startswith('# provenance\n'))
        self.assertEqual(csv.field_size_limit(), limit)

    def test_invalid_counts_and_columns_are_rejected_during_postprocessing(self):
        with tempfile.TemporaryDirectory() as temporary:
            source, output = Path(temporary) / 'source.txt', Path(temporary) / 'final.txt'
            header = '\t'.join((*ANNOTATION_COLUMNS, 'reads.sam')) + '\n'
            for row in ('gene\tchr1\t1\t20\t+\t20\t-1\n',
                        'gene\tchr1\t1\t20\t+\t20\t3\textra\n'):
                source.write_text(header + row)
                with self.assertRaises(ValueError):
                    add_zero_columns(source, output, ['reads.sam'], ['reads.sam'], ANNOTATION_COLUMNS)


@unittest.skipUnless(shutil.which('featureCounts'), 'featureCounts is required')
class PairedFeatureCountsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='paired counts ')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.gff = self.root / 'genes.gff'
        self.gff.write_text('chr\ttest\tCDS\t1\t20\t.\t+\t.\tlocus_tag=g1\nchr\ttest\tCDS\t15\t30\t.\t+\t.\tlocus_tag=g1\nchr\ttest\tCDS\t50\t90\t.\t-\t.\tlocus_tag=g2\n')
        self.output = self.root / 'counts.txt'
        self.empty_sam = self.root / 'empty.sam'
        self.empty_bam = self.root / 'empty.bam'
        header = {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chr', 'LN': 100}]}
        with pysam.AlignmentFile(str(self.empty_sam), 'w', header=header):
            pass
        with pysam.AlignmentFile(str(self.empty_bam), 'wb', header=header):
            pass
        self.pairs = self.root / 'pairs.sam'
        with pysam.AlignmentFile(str(self.pairs), 'w', header=header) as output:
            for flag, start, mate, length in ((83, 15, 5, -20), (163, 5, 15, 20)):
                read = pysam.AlignedSegment()
                read.query_name = 'reverse_strand_pair'
                read.flag = flag
                read.reference_id = read.next_reference_id = 0
                read.reference_start, read.next_reference_start = start, mate
                read.template_length = length
                read.mapping_quality = 42
                read.query_sequence = 'A' * 10
                read.cigarstring = '10M'
                read.query_qualities = pysam.qualitystring_to_array('I' * 10)
                output.write(read)

    def count(self, inputs):
        # Keep command output out of unittest logs while executing the real tool.
        actual_run = subprocess.run
        def run(command, **kwargs):
            return actual_run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, **kwargs)
        with patch('paired_featurecounts.subprocess.run', side_effect=run) as called:
            count_features(self.gff, 'CDS', self.output, 1, inputs)
        with self.output.open() as handle:
            matrix = list(csv.reader((line for line in handle if not line.startswith('#')), delimiter='\t'))
        with Path(str(self.output) + '.summary').open() as handle:
            summary = list(csv.reader(handle, delimiter='\t'))
        self.assertEqual(matrix[0][6:], list(map(str, inputs)))
        self.assertEqual(summary[0][1:], list(map(str, inputs)))
        self.assertEqual([row[5] for row in matrix[1:]], ['30', '41'])
        return matrix, summary, called.call_args.args[0]

    def test_all_empty_sam_and_bam_preserve_zero_columns_and_annotation_lengths(self):
        matrix, summary, command = self.count([self.empty_sam, self.empty_bam])
        self.assertEqual([row[6:] for row in matrix[1:]], [['0', '0'], ['0', '0']])
        self.assertTrue(all(row[1:] == ['0', '0'] for row in summary[1:]))
        self.assertNotIn('-p', command)
        self.assertEqual(command[command.index('-s') + 1], '2')

    def test_mixed_empty_nonempty_input_preserves_order_and_counts_pairs_once(self):
        matrix, summary, command = self.count([self.empty_sam, self.pairs, self.empty_bam])
        self.assertEqual([row[6:] for row in matrix[1:]], [['0', '1', '0'], ['0', '0', '0']])
        self.assertEqual(next(row for row in summary if row[0] == 'Assigned'), ['Assigned', '0', '1', '0'])
        for option in ('-p', '--countReadPairs', '-B', '-C'):
            self.assertIn(option, command)
        self.assertEqual(command[command.index('-s') + 1], '2')
        self.assertNotIn(str(self.empty_sam), command)
        self.assertNotIn(str(self.empty_bam), command)

    def test_nonempty_library_keeps_paired_reverse_stranded_counting(self):
        matrix, _, command = self.count([self.pairs])
        self.assertEqual([row[6:] for row in matrix[1:]], [['1'], ['0']])
        self.assertIn('-p', command)

    def test_unpaired_record_after_a_pair_is_rejected_before_running_counter(self):
        with self.pairs.open('a') as handle:
            handle.write('unpaired\t0\tchr\t5\t42\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n')
        with patch('paired_featurecounts.subprocess.run') as run, self.assertRaisesRegex(ValueError, 'paired-end'):
            count_features(self.gff, 'CDS', self.output, 1, [self.pairs])
        run.assert_not_called()
        self.assertFalse(self.output.exists())

    def test_real_featurecounts_with_many_exons_exceeds_old_csv_field_limit(self):
        self.gff.write_text(''.join(
            f'chr\ttest\tCDS\t{2*i+1}\t{2*i+1}\t.\t+\t.\tlocus_tag=large_gene\n'
            for i in range(30000)))
        actual_run = subprocess.run
        def run(command, **kwargs):
            return actual_run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, **kwargs)
        with patch('paired_featurecounts.subprocess.run', side_effect=run):
            count_features(self.gff, 'CDS', self.output, 1, [self.empty_sam, self.empty_bam])
        with open_table(self.output) as rows:
            header, row = list(rows)
        self.assertEqual(header[6:], [str(self.empty_sam), str(self.empty_bam)])
        self.assertGreater(max(map(len, row[:6])), 131072)
        self.assertEqual(row[5:], ['30000', '0', '0'])


if __name__ == '__main__':
    unittest.main()
