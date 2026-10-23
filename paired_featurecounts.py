#!/usr/bin/env python3
"""Run reverse-stranded paired featureCounts, including zero-record libraries.

featureCounts 2.1.1 rejects a header-only library when -p is set because no
alignment record exists from which to infer pairedness. Nonempty inputs always
use -p --countReadPairs -B -C -s 2. Empty libraries receive zero columns from the
same featureCounts annotation matrix. When every library is empty, featureCounts
loads the annotations using one verified zero-record input without -p; no reads
are processed in this annotation-only invocation. No artificial reads are added.
"""

import argparse
import csv
from pathlib import Path
import subprocess
import tempfile

import pysam


ANNOTATION_COLUMNS = ('Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length')


def contains_pairs(path):
    """Verify every alignment is paired; return False only for zero records."""
    any_records = False
    mode = 'rb' if str(path).endswith('.bam') else 'r'
    with pysam.AlignmentFile(str(path), mode) as alignments:
        for read in alignments.fetch(until_eof=True):
            if not read.is_paired:
                raise ValueError(f'Expected only paired-end alignments in {path}: {read.query_name}')
            any_records = True
    return any_records


def add_zero_columns(source, destination, inputs, counted, annotation_columns):
    """Preserve featureCounts annotation rows and insert empty-library columns."""
    positions = {name: index for index, name in enumerate(counted)}
    with source.open(newline='') as reader_handle, destination.open('w', newline='') as writer_handle:
        for line in reader_handle:
            if line.startswith('#'):
                writer_handle.write(line)
            else:
                header = next(csv.reader([line], delimiter='\t'))
                break
        else:
            raise ValueError(f'featureCounts produced no table header: {source}')
        if header[:len(annotation_columns)] != list(annotation_columns) or header[len(annotation_columns):] != counted:
            raise ValueError(f'Unexpected featureCounts output columns in {source}')
        writer = csv.writer(writer_handle, delimiter='\t', lineterminator='\n')
        writer.writerow((*annotation_columns, *inputs))
        for row in csv.reader(reader_handle, delimiter='\t'):
            if len(row) != len(annotation_columns) + len(counted):
                raise ValueError(f'Unexpected featureCounts output row in {source}')
            values = row[len(annotation_columns):]
            writer.writerow((*row[:len(annotation_columns)], *(values[positions[name]] if name in positions else '0' for name in inputs)))


def count_features(gff, feature, output, threads, inputs):
    if threads < 1:
        raise ValueError('threads must be a positive integer')
    if not inputs:
        raise ValueError('At least one paired-end alignment input is required')
    names = [str(path) for path in inputs]
    if len({Path(path).resolve() for path in names}) != len(names):
        raise ValueError('Duplicate alignment input files')
    nonempty = [name for name in names if contains_pairs(name)]
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.paired_featurecounts-', dir=output.parent) as temporary:
        work = Path(temporary)
        counted = nonempty or names[:1]
        result = work / 'counts.txt'
        command = ['featureCounts', '-T', str(threads), '-a', str(gff), '-s', '2', '-g', 'locus_tag', '-t', feature]
        if nonempty:
            command += ['-p', '--countReadPairs', '-B', '-C']
        command += ['-o', str(result), *counted]
        subprocess.run(command, check=True)
        # All-empty runs also get a zero column for the template input itself.
        # It contains zero alignments, so its annotation-only result must be zero.
        add_zero_columns(result, work / 'final.txt', names, counted, ANNOTATION_COLUMNS)
        add_zero_columns(Path(str(result) + '.summary'), work / 'final.summary', names, counted, ('Status',))
        (work / 'final.txt').replace(output)
        (work / 'final.summary').replace(Path(str(output) + '.summary'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gff', type=Path, required=True)
    parser.add_argument('--feature', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('inputs', nargs='+', type=Path)
    args = parser.parse_args()
    try:
        count_features(args.gff, args.feature, args.output, args.threads, args.inputs)
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()
