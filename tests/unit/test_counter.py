from mbarq.counter import BarcodeCounter
from collections import Counter
from pathlib import Path
import pandas as pd
import pickle
import subprocess
import shlex
import json
import random
from test_utils import assert_files_are_same

TESTDATA = "./tests/test_files"
EXPDATA = "./tests/expected_outcomes/counting"
OUTDIR = "./tests/tmp"


def test_extract_barcodes_rbseq(tn5_structure):
    """Check counting barcodes in fasta file"""
    r1 = f'{TESTDATA}/count_test.fasta'
    seq_data = BarcodeCounter(r1, tn5_structure, output_dir=OUTDIR)
    seq_data._extract_barcodes()
    with open(f"{EXPDATA}/test_extract_barcodes_rbseq.json", 'r') as fh:
        expected_counted = Counter(json.load(fh))
    assert(expected_counted == seq_data.barcode_counter)


def test_merge_similar(tn5_structure):
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
    seq_data = BarcodeCounter(r1, tn5_structure, output_dir=OUTDIR)
    seq_data.barcode_counter = barcodes
    seq_data._merge_similar()
    expected_counter = Counter({'AAGCCCAATAAACCACT': 105, 'CTGACTGGCCGAATAGG': 50,
                                'GATATAGGCAACGACAT': 43, 'GTGCGGCGACCCTTGCG': 22,
                                'ACAGTGACGCTTTCGCC': 11, 'GTTGCCTAAACCTATTT': 3})
    assert(expected_counter == seq_data.barcode_counter)
    assert(seq_data.merged is True)


def test_merge_similar_dnaid2023(tn5_structure, count_test_data_tn5, tmpdir):
    r1, mfile = count_test_data_tn5
    seq_data = BarcodeCounter(r1, tn5_structure, mapping_file=mfile, output_dir=tmpdir)
    seq_data._extract_barcodes()
    seq_data._merge_similar()
    with open(f"{EXPDATA}/test_merge_similar_dnaid2023.json", 'r') as fh:
        expected_counted = Counter(json.load(fh))
    assert(expected_counted == seq_data.barcode_counter)


def test_annotate_barcodes(tn5_structure, count_test_data_tn5, tmpdir):
    # todo use the latest map produced by mbarq
    r1, mfile = count_test_data_tn5
    seq_data = BarcodeCounter(r1, tn5_structure, mapping_file=mfile, output_dir=tmpdir)
    seq_data._extract_barcodes()
    seq_data._merge_similar()
    seq_data._annotate_barcodes()
    expected_counts = f"{EXPDATA}/test_annotate_barcodes_dnaid2023.csv"
    out_counts = tmpdir.join("test_annotate_barcodes_dnaid2023.csv")
    seq_data.annotated_cnts.to_csv(out_counts)
    assert_files_are_same(expected_counts, out_counts)


def test_count_barcodes(tn5_structure, count_test_data_tn5):
    # todo use new mapping file
    r1, mfile = count_test_data_tn5
    seq_data = BarcodeCounter(r1, tn5_structure, mapping_file=mfile, output_dir=OUTDIR)
    seq_data.count_barcodes()
    expected_counts = f"{EXPDATA}/test_annotate_barcodes_dnaid2023.csv"
    actual_counts = f"{OUTDIR}/dnaid2023_12_test_mbarq_counts.csv"
    assert_files_are_same(expected_counts, actual_counts)


def test_count_barcodes_wish(wish_structure, count_test_data_wish, tmpdir):
    fq, mfile = count_test_data_wish
    seq_data = BarcodeCounter(fq, wish_structure, mapping_file=mfile, output_dir=tmpdir)
    seq_data.count_barcodes()
    expected_counts = f"{EXPDATA}/LibraryA_pilot2_mbarq_counts.csv"
    actual_counts = tmpdir.join("LibraryA_pilot2_mbarq_counts.csv")
    assert_files_are_same(expected_counts, actual_counts)


#
