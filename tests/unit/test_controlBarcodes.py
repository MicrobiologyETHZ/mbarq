from mbarq.analysis import ControlBarcodes
from collections import Counter
from pathlib import Path
import pandas as pd
import pandera.errors
import pytest
import pickle
import subprocess
import shlex
import json
import random
from test_utils import assert_files_are_same

TESTDATA = "./tests/test_files/analysis"
EXPDATA = "./tests/expected_outcomes/analysis"
OUTDIR = "./tests/tmp"

WT_BCS = ['TACCCAGAGCACACTCA', 'ATCCGCGTCACCGAAAA', 'ACAGAGCTCGGGAGTCT',
                   'ACTACAAGACTGGTTAA', 'AGATGCATGACTAGCTA', 'AGAATGACCCGGAGGCT',
                   'AGGAAGGCGACGAAATC', 'AGTCATCGATGCTATAT', 'TAAGTCCGGGCTAAGTC',
                    'AACAACACGGTAAGCAA', 'TATAACACCCCCGATTC', 'CTACGACAGGGACTTAA',
                   'GTGTATAGCAGGAACCC', 'CCGACGACTGATTGTCC', 'TCTCACGCAGCGTTTCG']

def test__read_control_file_3cols():

    r1 = f'{TESTDATA}/controls_3col.csv'
    cntrlBC = ControlBarcodes(r1, OUTDIR)
    cntrlBC.read_control_file()
    assert(len(cntrlBC.wt_barcodes == len(WT_BCS)))
    assert(all(cntrlBC.wt_barcodes == WT_BCS))


def test__read_control_file_2cols():

    r1 = f'{TESTDATA}/controls_2col.csv'
    cntrlBC = ControlBarcodes(r1, OUTDIR)
    cntrlBC.read_control_file()
    assert(len(cntrlBC.wt_barcodes == len(WT_BCS)))
    assert(all(cntrlBC.wt_barcodes == WT_BCS))


def test__read_control_file_1cols():

    r1 = f'{TESTDATA}/controls_1col.csv'
    cntrlBC = ControlBarcodes(r1, OUTDIR)
    cntrlBC.read_control_file()
    assert(len(cntrlBC.wt_barcodes == len(WT_BCS)))
    assert(all(cntrlBC.wt_barcodes == WT_BCS))


def test__read_control_file_3cols_no_wt():

    r1 = f'{TESTDATA}/controls_3col_no_wt.csv'
    cntrlBC = ControlBarcodes(r1, OUTDIR)
    with pytest.raises(SystemExit):
        cntrlBC._read_control_file()
