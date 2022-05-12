import pandera.errors
from mbarq.analysis import CountDataSet
from collections import Counter
from pathlib import Path
import pandas as pd
import pickle
import subprocess
import shlex
import json
import random
from test_utils import assert_files_are_same
import pytest

TMPDIR = Path("/Users/ansintsova/git_repos/mbarq/tests/tmp")
EXPDATA = Path("/Users/ansintsova/git_repos/mbarq/tests/expected_outcomes/merging")


def test__merge_count_files(merge_test_data_tn5, tmpdir):
    count_files, gene_name = merge_test_data_tn5
    name = 'TestingMerge'
    cds = CountDataSet(count_files=count_files, name=name, gene_column_name=gene_name, output_dir=tmpdir)
    fdf = cds._merge_count_files()
    expected_csv = EXPDATA/"merge_tn5_no_validation.csv"
    out_file = tmpdir.join("merge_tn5_no_validation.csv")
    fdf.to_csv(out_file, index=False)
    assert_files_are_same(expected_csv, out_file)


def test__validate_count_table(merge_test_data_tn5, tmpdir):
    count_files, _ = merge_test_data_tn5
    name = 'TestingValidate'
    df1 = pd.DataFrame({
        "barcode": ['ACGT', 'ATTGC', 'ACCG', 'AAAA', 'CCC'],
        "gene": ["value_1", "value_2", "value_3", "value_2", "value_1"],
        "column2": [-1.3, -1.4, -2.9, -10.1, -20.4],
        "column3": [-1.3, -1.4, -2.9, -10.1, -20.4],
        "column4": [-1.3, -1.4, -2.9, 10, -20.4],

    })
    gene1 = "gene"
    cds = CountDataSet(count_files=count_files, name=name, gene_column_name=gene1, output_dir=tmpdir)
    _, sample_cols = cds._validate_count_table(df1)
    assert(Counter(sample_cols) == Counter(['column2', 'column3', 'column4']))

    df2 = pd.DataFrame({
        "barcode": ['ACGT', 'ATTGC', 'ACCG', 'AAAA', 'CCC'],
        "column2": [-1.3, -1.4, -2.9, -10.1, -20.4],
        "column3": [-1.3, -1.4, -2.9, -10.1, -20.4],
        "column4": [-1.3, -1.4, -2.9, 10, -20.4],

    })
    gene2 =''
    cds = CountDataSet(count_files=count_files, name=name, gene_column_name=gene2, output_dir=tmpdir)
    _, sample_cols = cds._validate_count_table(df2)
    assert (list(sample_cols) == ['column2', 'column3', 'column4'])

    with pytest.raises(SystemExit):
        df3 = pd.DataFrame({
            "barcode": ['ACGT', 'ATTGC', 'ACCG', 'AAAA', 'CCC'],
            "column2": [-1.3, -1.4, -2.9, -10.1, -20.4],
            "column3": [-1.3, -1.4, -2.9, -10.1, -20.4],
            "column4": [-1.3, -1.4, '', 10, -20.4],

        })
        gene3 = ''
        cds = CountDataSet(count_files=count_files, name=name, gene_column_name=gene3, output_dir=tmpdir)
        cds._validate_count_table(df3)


def test_create_count_table(merge_test_data_tn5, tmpdir):
    count_files, gene_name = merge_test_data_tn5
    name = 'TestCreateCountTable'
    cds = CountDataSet(count_files=count_files, name=name, gene_column_name=gene_name, output_dir=tmpdir)
    cds.create_count_table()
    expected_csv = EXPDATA / "merge_tn5_no_validation.csv"
    out_file = tmpdir.join("TestCreateCountTable_mbarq_merged_counts.csv")
    assert_files_are_same(expected_csv, out_file)

    count_files, gene_name = EXPDATA / "merge_tn5_no_validation.csv", 'Name'
    cds = CountDataSet(count_files=count_files, name=name, gene_column_name=gene_name, output_dir=tmpdir)
    cds.create_count_table()
    expected_csv = EXPDATA / "merge_tn5_no_validation.csv"
    out_file = tmpdir.join("TestCreateCountTable_mbarq_merged_counts.csv")
    assert_files_are_same(expected_csv, out_file)
