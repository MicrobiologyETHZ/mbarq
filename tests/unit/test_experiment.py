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

# Testing controlBarcodes

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


# Testing sampleData

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


def test_sample_data_fail(analysis_test_data_tn5, tmpdir):
    _,_,_,_,_,sample_file,_ = analysis_test_data_tn5
    # batch column does not exist
    sd = SampleData(sample_file, 'Test', 'day', 'd0', 'batch', tmpdir)
    with pytest.raises(SystemExit):
        sd.read_sample_data_csv()
    # baseline does not exist
    sd = SampleData(sample_file, 'Test', 'day', 'd9', 'experiment', tmpdir)
    with pytest.raises(SystemExit):
        sd.read_sample_data_csv()
    sd = SampleData(sample_file, 'Test', 'cucumber', 'd0', 'experiment', tmpdir)
    with pytest.raises(SystemExit):
        sd.read_sample_data_csv()


# Testing experiment

def test_get_good_samples(analysis_test_data_tn5, tmpdir):
    _,_,_,controls,count_file,sample_file,_ = analysis_test_data_tn5
    name = 'TestExp'
    # cds = CountDataSet(count_file, name=name, gene_column_name='Name', output_dir=tmpdir)
    # cds.create_count_table()
    # cds.calculate_cpms()
    # cbars = ControlBarcodes(controls, output_dir = tmpdir)
    # cbars.read_control_file()
    # sd = SampleData(sample_file, name,'day', 'd0', 'experiment', tmpdir )
    # sd.read_sample_data_csv()
    exp = Experiment(count_file, sample_file, controls, name, 'Name', 'day', 'd0', 'experiment', 0.8, tmpdir)
    exp._get_good_samples()
    expected_good_samples = ['dnaid1315_10', 'dnaid1315_107', 'dnaid1315_117', 'dnaid1315_124', 'dnaid1315_128',
                             'dnaid1315_129', 'dnaid1315_131', 'dnaid1315_136', 'dnaid1315_17', 'dnaid1315_18',
                             'dnaid1315_19', 'dnaid1315_20', 'dnaid1315_28', 'dnaid1315_40', 'dnaid1315_42',
                             'dnaid1315_50', 'dnaid1315_52', 'dnaid1315_66', 'dnaid1315_81', 'dnaid1315_90',
                             'dnaid1315_92', 'dnaid1315_94', 'dnaid1315_96']
    assert exp.good_samples == expected_good_samples

def test_prepare_mageck_dataset(analysis_test_data_tn5, tmpdir, dnaid1315_expected_outcomes):
    _, _, _, controls, count_file, sample_file, _ = analysis_test_data_tn5
    name = "test_prepare_mageck_dataset"
    exp = Experiment(count_file, sample_file, controls, name, 'Name', 'day', 'd0', 'experiment', 0.8, tmpdir)
    exp._get_good_samples()
    exp.prepare_mageck_dataset()
    expected_counts = dnaid1315_expected_outcomes /"test_prepare_mageck_dataset_count.txt"
    expected_sd = dnaid1315_expected_outcomes/"test_prepare_mageck_dataset_batch.txt"
    assert exp.batch_file.is_file()
    assert exp.count_file.is_file()
    assert_files_are_same(expected_counts, exp.count_file)
    assert_files_are_same(expected_sd, exp.batch_file)


def test_batch_correction(analysis_test_data_tn5, tmpdir, dnaid1315_expected_outcomes):
    _, _, _, controls, count_file, sample_file, _ = analysis_test_data_tn5
    name = "test_batch_correction"
    exp = Experiment(count_file, sample_file, controls, name, 'Name', 'day', 'd0', 'experiment', 0.8, tmpdir)
    exp._get_good_samples()
    exp.prepare_mageck_dataset()
    exp.batch_correct()
    expected_count = dnaid1315_expected_outcomes/"test_batch_correction_count.batchcorrected.txt"
    assert_files_are_same(exp.count_file, expected_count)


def test_get_contrast_samples(analysis_test_data_tn5, tmpdir, dnaid1315_expected_outcomes):
    _, _, _, controls, count_file, sample_file, _ = analysis_test_data_tn5
    name = "test_get_contrast_samples"
    exp = Experiment(count_file, sample_file, controls, name, 'Name', 'day', 'd0', 'experiment', 0.8, tmpdir)
    exp._get_good_samples()
    treatment = 'd1'
    controls, treats = exp.get_contrast_samples(treatment)
    assert controls == 'dnaid1315_10,dnaid1315_81,dnaid1315_107'
    assert treats == 'dnaid1315_17,dnaid1315_18,dnaid1315_19,dnaid1315_90,dnaid1315_117,dnaid1315_124'


def test_run_mageck(analysis_test_data_tn5, tmpdir, dnaid1315_expected_outcomes):
    _, _, _, controls, count_file, sample_file, _ = analysis_test_data_tn5
    name = "test_run_mageck"
    exp = Experiment(count_file, sample_file, controls, name, 'Name', 'day', 'd0', 'experiment', 0.8, tmpdir)
    exp._get_good_samples()
    exp.count_file = dnaid1315_expected_outcomes/"test_batch_correction_count.batchcorrected.txt"
    exp.batch_file = dnaid1315_expected_outcomes/"test_prepare_mageck_dataset_batch.txt"
    exp.write_control_barcodes_to_file()
    treatment = 'd1'
    controls, treats = exp.get_contrast_samples(treatment)
    run_name = "test_run_mageck_d1_vs_d0"
    exp.run_mageck(treats, controls, run_name)
    expected_gene_summary = dnaid1315_expected_outcomes/"test_run_mageck.gene_summary.txt"
    actual_gene_summary = tmpdir.join(f"{run_name}.gene_summary.txt")
    assert_files_are_same(actual_gene_summary, expected_gene_summary)

def test_run_all_contrasts(analysis_test_data_tn5, tmpdir, dnaid1315_expected_outcomes):
    _, _, _, controls, count_file, sample_file, _ = analysis_test_data_tn5
    name = "test_run_all_contrasts"
    exp = Experiment(count_file, sample_file, controls, name, 'Name', 'day', 'd0', 'experiment', 0.8, tmpdir)
    exp._get_good_samples()
    exp.count_file = dnaid1315_expected_outcomes / "test_batch_correction_count.batchcorrected.txt"
    exp.batch_file = dnaid1315_expected_outcomes / "test_prepare_mageck_dataset_batch.txt"
    exp.write_control_barcodes_to_file()
    exp.run_all_contrasts()
    expected_gene_summary = dnaid1315_expected_outcomes / "test_run_all_contrasts_d1_vs_d0.gene_summary.txt"
    actual_gene_summary = tmpdir.join("test_run_all_contrasts_d1_vs_d0.gene_summary.txt")
    assert_files_are_same(actual_gene_summary, expected_gene_summary)
    expected_gene_summary = dnaid1315_expected_outcomes / "test_run_all_contrasts_d2_vs_d0.gene_summary.txt"
    actual_gene_summary = tmpdir.join("test_run_all_contrasts_d2_vs_d0.gene_summary.txt")
    assert_files_are_same(actual_gene_summary, expected_gene_summary)
    expected_gene_summary = dnaid1315_expected_outcomes / "test_run_all_contrasts_d3_vs_d0.gene_summary.txt"
    actual_gene_summary = tmpdir.join("test_run_all_contrasts_d3_vs_d0.gene_summary.txt")
    assert_files_are_same(actual_gene_summary, expected_gene_summary)
    expected_gene_summary = dnaid1315_expected_outcomes / "test_run_all_contrasts_d4_vs_d0.gene_summary.txt"
    actual_gene_summary = tmpdir.join("test_run_all_contrasts_d4_vs_d0.gene_summary.txt")
    assert_files_are_same(actual_gene_summary, expected_gene_summary)


def test_process_results(analysis_test_data_tn5, tmpdir, dnaid1315_expected_outcomes):
    _, _, _, controls, count_file, sample_file, _ = analysis_test_data_tn5
    name = "test_process_results"
    exp = Experiment(count_file, sample_file, controls, name, 'Name', 'day', 'd0', 'experiment', 0.8, tmpdir)
    exp._get_good_samples()
    exp.count_file = dnaid1315_expected_outcomes / "test_batch_correction_count.batchcorrected.txt"
    exp.batch_file = dnaid1315_expected_outcomes / "test_prepare_mageck_dataset_batch.txt"
    exp.write_control_barcodes_to_file()
    exp.run_all_contrasts()
    exp.process_results()
    expected_rra = dnaid1315_expected_outcomes/"test_process_results_rra_results_no_index.csv"
    actual_rra = tmpdir.join("test_process_results_rra_results.csv")
    assert_files_are_same(actual_rra, expected_rra)


def test_run_experiment(analysis_test_data_tn5, tmpdir, dnaid1315_expected_outcomes):
    _, _, _, controls, count_file, sample_file, _ = analysis_test_data_tn5
    name = "test_run_experiment"
    exp = Experiment(count_file, sample_file, controls, name, 'Name', 'day', 'd0', 'experiment', 0.8, tmpdir)
    exp.run_experiment()
    expected_rra = dnaid1315_expected_outcomes / "test_process_results_rra_results_no_index.csv"
    actual_rra = tmpdir.join("test_run_experiment_rra_results.csv")
    assert_files_are_same(actual_rra, expected_rra)

OUTDIR= "/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/ansintsova/Projects_NCCR/hardt/nguyenb/tnseq/scratch/tmp"

def test_run_experiment_no_batch(analysis_test_data_tn5, tmpdir, dnaid1315_expected_outcomes):
    _, _, _, controls, count_file, sample_file, _ = analysis_test_data_tn5
    name = "test_run_experiment"
    exp = Experiment(count_file, sample_file, controls, name, 'Name', 'day', 'd0', '', 0.8, tmpdir)
    exp.run_experiment()
    expected_rra = dnaid1315_expected_outcomes / "test_run_experiment_rra_results_no_batch_no_index.csv"
    actual_rra = tmpdir.join("test_run_experiment_rra_results.csv")
    assert_files_are_same(actual_rra, expected_rra)


def test_run_experiment_no_control(analysis_test_data_tn5, tmpdir, dnaid1315_expected_outcomes):
    _, _, _, controls, count_file, sample_file, _ = analysis_test_data_tn5
    name = "test_run_experiment"
    exp = Experiment(count_file, sample_file, '', name, 'Name', 'day', 'd0', 'experiment', 0.8, tmpdir)
    exp.run_experiment()
    expected_rra = dnaid1315_expected_outcomes / "test_run_experiment_rra_results_no_control_no_index.csv"
    actual_rra = tmpdir.join("test_run_experiment_rra_results.csv")
    assert_files_are_same(actual_rra, expected_rra)
