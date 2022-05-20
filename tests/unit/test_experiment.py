from mbarq.analysis import CountDataSet, ControlBarcodes, Experiment, SampleData
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


def test_sample_data():
    sample_file = f"{TESTDATA}/example_sample_data.csv"
    sd = SampleData(sample_file, 'Test', 'da', 'd1', 'experiment', OUTDIR)
    sd.read_sample_data_csv()
    print(sd.sampleData)

def test_calculate_cpms():
    count_file = f"{TESTDATA}/example_mbarq_merged_counts.csv"
    countDS = CountDataSet(count_file, name='test', gene_column_name='Name', output_dir=OUTDIR)
    countDS.create_count_table()
    countDS.calculate_cpms()
    controlBCs = ControlBarcodes(f"{TESTDATA}/controls_3col.csv", output_dir = OUTDIR)
    controlBCs.read_control_file()
    exp = Experiment(countDS, controlBCs)
    exp._get_good_samples()
    print(exp.good_samples)


def test_calculate_cpms_2col():
    count_file = f"{TESTDATA}/example_mbarq_merged_counts.csv"
    countDS = CountDataSet(count_file, name='test', gene_column_name='Name', output_dir=OUTDIR)
    countDS.create_count_table()
    countDS.calculate_cpms()
    controlBCs = ControlBarcodes(f"{TESTDATA}/controls_2col_short.csv", output_dir = OUTDIR)
    controlBCs.read_control_file()
    sample_file = f"{TESTDATA}/example_sample_data.csv"
    sd = SampleData(sample_file, 'Test', 'day', 'd1', 'experiment', OUTDIR)
    sd.read_sample_data_csv()
    exp = Experiment('Test', countDS, controlBCs, sd)
    exp._validate_experiment()
    exp._get_good_samples()
    exp.prepare_mageck_dataset()
    print(exp.good_samples)


def test_batch_correction():
    count_file = f"{TESTDATA}/example_mbarq_merged_counts.csv"
    countDS = CountDataSet(count_file, name='test', gene_column_name='Name', output_dir=OUTDIR)
    countDS.create_count_table()
    countDS.calculate_cpms()
    controlBCs = ControlBarcodes(f"{TESTDATA}/controls_2col_short.csv", output_dir=OUTDIR)
    controlBCs.read_control_file()
    sample_file = f"{TESTDATA}/example_sample_data.csv"
    sd = SampleData(sample_file, 'Test', 'day', 'd1', 'experiment', OUTDIR)
    sd.read_sample_data_csv()
    exp = Experiment('Test', countDS, controlBCs, sd, output_dir=OUTDIR)
    print(exp.name)
    exp._validate_experiment()
    exp._get_good_samples()
    exp.prepare_mageck_dataset()
    exp.batch_correct()

#test_calculate_cpms_2col()
#test_sample_data()
test_batch_correction()
