from mbarq.counter import BarcodeCounter
from collections import Counter
from pathlib import Path
import pandas as pd
import pickle
import subprocess
import shlex
import json
import random


TESTDATA = "./tests/test_files"
EXPDATA = "./tests/expected_outcomes/counting"
OUTDIR = "./tests/tmp"



#  todo Add this code to the conftest.py?


def get_structure(experiment='rbseq'):
    if experiment == 'rbseq':
        #return 'GTGTATAAGAGACAG:17:13:before'
        return 'B17N13GTGTATAAGAGACAG'
    elif experiment == 'wish':
        return 'GGAGGTTCACAATGTGGGAGGTCAB40'
    else:
        return None


def capture(command_str):
    command = shlex.split(command_str)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,)
    out, err = proc.communicate()
    return out, err, proc


def assert_files_are_same(file1, file2, verbose=False):
    cmd_str = f'cmp {file1} {file2}'
    out, err, proc = capture(cmd_str)
    if verbose:
        print(out)
        print(err)
    assert proc.returncode == 0


def get_test_file(experiment='rbseq'):
    if experiment == 'rbseq':
        return f'{TESTDATA}/dnaid2023_12_test.fasta', f'{TESTDATA}/ref/library_13_1.barcode_map.annotated.csv'
    elif experiment == 'wish':
        return f'{TESTDATA}/LibraryA_pilot2.fq.gz', f'{TESTDATA}/ref/20210520_BarcodeList.csv'
    else:
        return None


def test_extract_barcodes_rbseq():
    """Check counting barcodes in fasta file"""
    r1 = f'{TESTDATA}/count_test.fasta'
    structure = get_structure()
    seq_data = BarcodeCounter(r1, structure, output_dir=OUTDIR)
    seq_data._extract_barcodes()
    with open(f"{EXPDATA}/test_extract_barcodes_rbseq.json", 'r') as fh:
        expected_counted = Counter(json.load(fh))
    assert(expected_counted == seq_data.barcode_counter)


def test_merge_similar():
    """

       Final resutls:

       AAGCCCAATAAACCACT 105
       CTGACTGGCCGAATAGG 50
       GATATAGGCAACGACAT 43
       GTGCGGCGACCCTTGCG 22
       ACAGTGACGCTTTCGCC 11
       GTTGCCTAAACCTATTT 3

       :return:
       """
    #random.seed(42)
    #barcodes = ["".join([random.choice(['A', 'C', 'G', 'T']) for _ in range(17)]) for _ in range(6)]

    barcodes = Counter(
        {'AAGCCCAATAAACCACT': 100, 'CTGACTGGCCGAATAGG': 50, 'GATATAGGCAACGACAT': 40,
         'GTGCGGCGACCCTTGCG': 20, 'ACAGTGACGCTTTCGCC': 10, 'GTTGCCTAAACCTATTT': 3,
         'TTGCCCAATAAACCACT': 5, 'GATATCGGCAACGACAT': 3, 'GTGCGCCGAACCTTGCG': 2,
         'ACAGTGAGGCTTTCGCC': 1})

    r1 = f'{TESTDATA}/count_test.fasta'
    structure = get_structure()
    seq_data = BarcodeCounter(r1, structure, output_dir=OUTDIR)
    seq_data.barcode_counter = barcodes
    seq_data._merge_similar()
    expected_counter = Counter({'AAGCCCAATAAACCACT': 105, 'CTGACTGGCCGAATAGG': 50,
                                'GATATAGGCAACGACAT': 43, 'GTGCGGCGACCCTTGCG': 22,
                                'ACAGTGACGCTTTCGCC': 11, 'GTTGCCTAAACCTATTT': 3})
    assert(expected_counter == seq_data.barcode_counter)
    assert(seq_data.merged is True)


def test_merge_similar_dnaid2023():
    r1, mfile = get_test_file(experiment='rbseq')
    structure = get_structure()
    seq_data = BarcodeCounter(r1, structure, mapping_file=mfile, output_dir=OUTDIR)
    seq_data._extract_barcodes()
    seq_data._merge_similar()
    with open(f"{EXPDATA}/test_merge_similar_dnaid2023.json", 'r') as fh:
        expected_counted = Counter(json.load(fh))
    assert(expected_counted == seq_data.barcode_counter)


def test_annotate_barcodes():
    r1, mfile = get_test_file(experiment='rbseq')
    structure = get_structure()
    seq_data = BarcodeCounter(r1, structure, mapping_file=mfile, output_dir=OUTDIR)
    seq_data._extract_barcodes()
    seq_data._merge_similar()
    seq_data._annotate_barcodes()
    expected_counts = f"{EXPDATA}/test_annotate_barcodes_dnaid2023.csv"
    seq_data.annotated_cnts.to_csv(f"{OUTDIR}/test_annotate_barcodes_dnaid2023.csv")
    assert_files_are_same(expected_counts, f"{OUTDIR}/test_annotate_barcodes_dnaid2023.csv")


def test_count_barcodes():
    r1, mfile = get_test_file(experiment='rbseq')
    structure = get_structure()
    seq_data = BarcodeCounter(r1, structure, mapping_file=mfile, output_dir=OUTDIR)
    seq_data.count_barcodes()
    expected_counts = f"{EXPDATA}/test_annotate_barcodes_dnaid2023.csv"
    actual_counts = f"{OUTDIR}/dnaid2023_12_test_mbarq_counts.csv"
    assert_files_are_same(expected_counts, actual_counts)

def test_count_barcodes_wish():
    fq, mfile = get_test_file(experiment='wish')
    structure = get_structure(experiment='wish')
    seq_data = BarcodeCounter(fq, structure, mapping_file=mfile, output_dir=OUTDIR)
    seq_data.count_barcodes()
    expected_counts = f"{EXPDATA}/LibraryA_pilot2_mbarq_counts.csv"
    actual_counts = f"{OUTDIR}/LibraryA_pilot2_mbarq_counts.csv"
    assert_files_are_same(expected_counts, actual_counts)



if __name__ == '__main__':
    test_merge_similar_dnaid2023()