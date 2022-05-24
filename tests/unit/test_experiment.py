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

# TESTDATA = "./tests/test_files/analysis"
# EXPDATA = "./tests/expected_outcomes/analysis"
# OUTDIR = "./tests/tmp"

def test_sample_data(analysis_test_data_tn5, tmpdir):
    _,_,_,_,_,sample_file,_ = analysis_test_data_tn5
    sd = SampleData(sample_file, 'Test', 'day', 'd0', 'experiment', tmpdir)
    sd.read_sample_data_csv()
    expected_sampleIDs = ['dnaid1315_10', 'dnaid1315_17', 'dnaid1315_18', 'dnaid1315_19',
                          'dnaid1315_20', 'dnaid1315_28', 'dnaid1315_40', 'dnaid1315_42',
                          'dnaid1315_50', 'dnaid1315_52', 'dnaid1315_53', 'dnaid1315_63',
                          'dnaid1315_66', 'dnaid1315_81', 'dnaid1315_90', 'dnaid1315_92',
                          'dnaid1315_94', 'dnaid1315_96', 'dnaid1315_107', 'dnaid1315_117',
                          'dnaid1315_124', 'dnaid1315_128', 'dnaid1315_129', 'dnaid1315_131',
                          'dnaid1315_136', 'dnaid1315_137', 'dnaid1315_139']
    expected_contrasts = ['d1', 'd2', 'd3', 'd4']
    assert (sd.contrasts == expected_contrasts)
    assert (sd.sampleIDs == expected_sampleIDs)

def test_get_good_samples(analysis_test_data_tn5, tmpdir):
    _,_,_,controls,count_file,sample_file,_ = analysis_test_data_tn5
    name = 'TestExp'
    cds = CountDataSet(count_file, name=name, gene_column_name='Name', output_dir=tmpdir)
    cds.create_count_table()
    cds.calculate_cpms()
    cbars = ControlBarcodes(controls, output_dir = tmpdir)
    cbars.read_control_file()
    sd = SampleData(sample_file, name,'day', 'd0', 'experiment', tmpdir )
    sd.read_sample_data_csv()
    exp = Experiment(name, cds, cbars, sd, output_dir=tmpdir)
    exp._get_good_samples()
    expected_good_samples = ['dnaid1315_10', 'dnaid1315_107', 'dnaid1315_117', 'dnaid1315_124', 'dnaid1315_128',
                             'dnaid1315_129', 'dnaid1315_131', 'dnaid1315_136', 'dnaid1315_17', 'dnaid1315_18',
                             'dnaid1315_19', 'dnaid1315_20', 'dnaid1315_28', 'dnaid1315_40', 'dnaid1315_42',
                             'dnaid1315_50', 'dnaid1315_52', 'dnaid1315_66', 'dnaid1315_81', 'dnaid1315_90',
                             'dnaid1315_92', 'dnaid1315_94', 'dnaid1315_96']
    assert exp.good_samples == expected_good_samples


def test_prepare_mageck_dataset(analysis_test_data_tn5, tmpdir):
    _, _, _, controls, count_file, sample_file, _ = analysis_test_data_tn5
    name = "test_calculate_cpms_2col"
    cds = CountDataSet(count_file, name=name, gene_column_name='Name', output_dir=tmpdir)
    cds.create_count_table()
    cds.calculate_cpms()
    cbars = ControlBarcodes(controls, output_dir = tmpdir)
    cbars.read_control_file()
    sd = SampleData(sample_file, 'Test', 'day', 'd0', 'experiment', tmpdir)
    sd.read_sample_data_csv()
    exp = Experiment(name, cds, cbars, sd, output_dir=tmpdir)
    exp._get_good_samples()
    exp.prepare_mageck_dataset()
    assert exp.batch_file.is_file()
    assert exp.count_file.is_file()



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
# test_batch_correction()
