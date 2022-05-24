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

WT_BCS = ['TACCCAGAGCACACTCA', 'ATCCGCGTCACCGAAAA', 'ACAGAGCTCGGGAGTCT',
                   'ACTACAAGACTGGTTAA', 'AGATGCATGACTAGCTA', 'AGAATGACCCGGAGGCT',
                   'AGGAAGGCGACGAAATC', 'AGTCATCGATGCTATAT', 'TAAGTCCGGGCTAAGTC',
                    'AACAACACGGTAAGCAA', 'TATAACACCCCCGATTC', 'CTACGACAGGGACTTAA',
                   'GTGTATAGCAGGAACCC', 'CCGACGACTGATTGTCC', 'TCTCACGCAGCGTTTCG']

def test__read_control_file_3cols(analysis_test_data_tn5, tmpdir):
    _,_,r1,_,_,_,_ = analysis_test_data_tn5
    cntrlBC = ControlBarcodes(r1, tmpdir)
    cntrlBC.read_control_file()
    assert(len(cntrlBC.wt_barcodes == len(WT_BCS)))
    assert(all(cntrlBC.wt_barcodes == WT_BCS))


def test__read_control_file_2cols(analysis_test_data_tn5, tmpdir):
    _,r1,_,_,_,_,_ = analysis_test_data_tn5
    cntrlBC = ControlBarcodes(r1, tmpdir)
    cntrlBC.read_control_file()
    assert(len(cntrlBC.wt_barcodes == len(WT_BCS)))
    assert(all(cntrlBC.wt_barcodes == WT_BCS))


def test__read_control_file_1cols(analysis_test_data_tn5, tmpdir):
    r1,_,_,_,_,_,_ = analysis_test_data_tn5
    cntrlBC = ControlBarcodes(r1, tmpdir)
    cntrlBC.read_control_file()
    assert(len(cntrlBC.wt_barcodes == len(WT_BCS)))
    assert(all(cntrlBC.wt_barcodes == WT_BCS))


def test__read_control_file_3cols_no_wt(analysis_test_data_tn5, tmpdir):
    _,_,_,_,_,_,r1 = analysis_test_data_tn5
    cntrlBC = ControlBarcodes(r1, tmpdir)
    with pytest.raises(SystemExit):
        cntrlBC.read_control_file()
