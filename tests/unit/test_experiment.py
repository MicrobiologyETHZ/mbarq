from mbarq.analysis import CountDataSet, ControlBarcodes, Experiment
from collections import Counter
from pathlib import Path
import pandas as pd
import pretty_errors
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


def test_calculate_cpms():
    count_file = f"{TESTDATA}/example_mbarq_merged_counts.csv"
    countDS = CountDataSet(count_file, name='test', gene_column_name='Name', output_dir=OUTDIR)
    countDS.create_count_table()
    countDS.calculate_cpms()
    controlBCs = ControlBarcodes(f"{TESTDATA}/controls_3col.csv", output_dir = OUTDIR)
    controlBCs.read_control_file()
    exp = Experiment(countDS, controlBCs)
    exp._get_wt_bc_counts()
    print(exp.wt_bc_conc_counts.head())


test_calculate_cpms()