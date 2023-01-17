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

# TESTDATA = "./tests/test_files"
# EXPDATA = "./tests/expected_outcomes/counting"
# OUTDIR = "./tests/tmp"


def test_extract_barcodes_rbseq(count_test_data_tn5, tn5_structure, tmpdir, dnaid1315_expected_outcomes):
    """Check counting barcodes in fasta file"""
    r1, _, _ = count_test_data_tn5
    seq_data = BarcodeCounter(r1, tn5_structure, output_dir=tmpdir)
    seq_data._extract_barcodes()
    with open(dnaid1315_expected_outcomes/"test_extract_barcodes_rbseq.json", 'r') as fh:
        expected_counted = Counter(json.load(fh))
    assert(expected_counted == seq_data.barcode_counter)


def test_merge_similar(count_test_data_tn5, tn5_structure, tmpdir):
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

    r1, _, _ = count_test_data_tn5
    seq_data = BarcodeCounter(r1, tn5_structure, output_dir=tmpdir)
    seq_data.barcode_counter = barcodes
    seq_data._merge_similar()
    expected_counter = Counter({'AAGCCCAATAAACCACT': 105, 'CTGACTGGCCGAATAGG': 50,
                                'GATATAGGCAACGACAT': 43, 'GTGCGGCGACCCTTGCG': 22,
                                'ACAGTGACGCTTTCGCC': 11, 'GTTGCCTAAACCTATTT': 3})
    assert(expected_counter == seq_data.barcode_counter)
    assert(seq_data.merged is True)


def test_merge_similar_dnaid1315(tn5_structure, count_test_data_tn5, tmpdir, dnaid1315_expected_outcomes):
    _, r1, mfile = count_test_data_tn5
    seq_data = BarcodeCounter(r1, tn5_structure, mapping_file=mfile, output_dir=tmpdir)
    seq_data._extract_barcodes()
    seq_data._merge_similar()
    with open(dnaid1315_expected_outcomes/"test_merge_similar_dnaid1315_124.json", 'r') as fh:
        expected_counted = Counter(json.load(fh))
    assert(expected_counted == seq_data.barcode_counter)


def test_annotate_barcodes(tn5_structure, count_test_data_tn5, tmpdir, dnaid1315_expected_outcomes):
    # todo use the latest map produced by mbarq
    _, r1, mfile = count_test_data_tn5
    seq_data = BarcodeCounter(r1, tn5_structure, mapping_file=mfile, output_dir=tmpdir)
    seq_data._extract_barcodes()
    seq_data._merge_similar()
    seq_data._annotate_barcodes(filter_low=True)
    expected_counts = dnaid1315_expected_outcomes/"test_annotate_barcodes_dnaid1315_124.csv"
    out_counts = tmpdir.join("test_annotate_barcodes_dnaid1315_124.csv")
    seq_data.annotated_cnts.to_csv(out_counts)
    assert_files_are_same(expected_counts, out_counts)

# todo same without filtering low counts

def test_count_barcodes(tn5_structure, count_test_data_tn5, dnaid1315_expected_outcomes, tmpdir):
    # todo use new mapping file
    _, r1, mfile = count_test_data_tn5
    seq_data = BarcodeCounter(r1, tn5_structure, mapping_file=mfile, output_dir= tmpdir)
    seq_data.count_barcodes(filter_low=True)
    expected_counts = dnaid1315_expected_outcomes/"dnaid1315_124_subsample_mbarq_counts.csv"
    actual_counts = tmpdir.join("dnaid1315_124_subsample_mbarq_counts.csv")
    assert_files_are_same(expected_counts, actual_counts)


# def test_count_barcodes_wish(wish_structure, count_test_data_wish, tmpdir):
#     fq, mfile = count_test_data_wish
#     seq_data = BarcodeCounter(fq, wish_structure, mapping_file=mfile, output_dir=tmpdir)
#     seq_data.count_barcodes()
#     expected_counts = f"{EXPDATA}/LibraryA_pilot2_mbarq_counts.csv"
#     actual_counts = tmpdir.join("LibraryA_pilot2_mbarq_counts.csv")
#     assert_files_are_same(expected_counts, actual_counts)


#
