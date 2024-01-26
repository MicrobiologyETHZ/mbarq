from mbarq.analysis import CountDataSet
from collections import Counter
from pathlib import Path
import pandas as pd
import pytest


def test__merge_count_files(merge_test_data, tmpdir):
    count_files, gene_name = merge_test_data
    name = "TestingMerge"
    cds = CountDataSet(
        count_files=count_files,
        name=name,
        gene_column_name=gene_name,
        output_dir=tmpdir,
    )
    fdf = cds._merge_count_files()
    expected_df = pd.DataFrame(
        {
            "barcode": {
                0: "AAAACCTCCCTGCCCAT",
                1: "ACGCAGACCCTCACTTT",
                2: "ACGCGCCAGACTTACGC",
                3: "AGAAGGAGAGCGAATAT",
                4: "AGATAACGAAACCACAC",
                5: "AGATCGCTGCTCGGGCG",
                6: "CTTCACTGTCATACGAA",
                7: "GAGCTAACCGATAACGG",
                8: "GATAGCTTGATGACGCA",
                9: "GTGACTCCGTCCAACAG",
                10: "TAAACAATGTACATAGA",
                11: "TATCGAACCACATCATA",
                12: "TCTGCAACGAGTTCAAC",
                13: "TGTTTTGGTAACGCTGC",
            },
            "Name": {
                0: "strB",
                1: "mobA",
                2: "traY",
                3: "hilC",
                4: "yjdB",
                5: "SL1344_4046",
                6: "garD",
                7: "SL1344_3745",
                8: "SL1344_3959",
                9: "hilD",
                10: "hilA",
                11: "hilD",
                12: "yhjV",
                13: "hilD",
            },
            "cf1": {
                0: 12571.0,
                1: 0.0,
                2: 19424.0,
                3: 0.0,
                4: 12611.0,
                5: 13661.0,
                6: 14052.0,
                7: 0.0,
                8: 14535.0,
                9: 0.0,
                10: 0.0,
                11: 18368.0,
                12: 15462.0,
                13: 26486.0,
            },
            "cf2": {
                0: 0.0,
                1: 13683.0,
                2: 18312.0,
                3: 0.0,
                4: 11722.0,
                5: 11964.0,
                6: 12879.0,
                7: 0.0,
                8: 12579.0,
                9: 0.0,
                10: 0.0,
                11: 17041.0,
                12: 15061.0,
                13: 23119.0,
            },
            "cf3": {
                0: 0.0,
                1: 0.0,
                2: 0.0,
                3: 21138.0,
                4: 0.0,
                5: 0.0,
                6: 20477.0,
                7: 17262.0,
                8: 18396.0,
                9: 28086.0,
                10: 19519.0,
                11: 51787.0,
                12: 22450.0,
                13: 59248.0,
            },
        }
    )
    assert fdf.equals(expected_df)


def test__validate_count_table(merge_test_data, tmpdir):
    count_files, _ = merge_test_data
    name = "TestingValidate"
    df1 = pd.DataFrame(
        {
            "barcode": ["ACGT", "ATTGC", "ACCG", "AAAA", "CCC"],
            "gene": ["value_1", "value_2", "value_3", "value_2", "value_1"],
            "column2": [-1.3, -1.4, -2.9, -10.1, -20.4],
            "column3": [-1.3, -1.4, -2.9, -10.1, -20.4],
            "column4": [-1.3, -1.4, -2.9, 10, -20.4],
        }
    )
    gene1 = "gene"
    cds = CountDataSet(
        count_files=count_files, name=name, gene_column_name=gene1, output_dir=tmpdir
    )
    _, sample_cols = cds._validate_count_table(df1)
    assert Counter(sample_cols) == Counter(["column2", "column3", "column4"])

    df2 = pd.DataFrame(
        {
            "barcode": ["ACGT", "ATTGC", "ACCG", "AAAA", "CCC"],
            "column2": [-1.3, -1.4, -2.9, -10.1, -20.4],
            "column3": [-1.3, -1.4, -2.9, -10.1, -20.4],
            "column4": [-1.3, -1.4, -2.9, 10, -20.4],
        }
    )
    gene2 = ""
    cds = CountDataSet(
        count_files=count_files, name=name, gene_column_name=gene2, output_dir=tmpdir
    )
    _, sample_cols = cds._validate_count_table(df2)
    assert list(sample_cols) == ["column2", "column3", "column4"]

    with pytest.raises(SystemExit):
        df3 = pd.DataFrame(
            {
                "barcode": ["ACGT", "ATTGC", "ACCG", "AAAA", "CCC"],
                "column2": [-1.3, -1.4, -2.9, -10.1, -20.4],
                "column3": [-1.3, -1.4, -2.9, -10.1, -20.4],
                "column4": [-1.3, -1.4, "", 10, -20.4],
            }
        )
        gene3 = ""
        cds = CountDataSet(
            count_files=count_files,
            name=name,
            gene_column_name=gene3,
            output_dir=tmpdir,
        )
        cds._validate_count_table(df3)


def test_create_count_table(merge_test_data, tmpdir):
    count_files, gene_name = merge_test_data
    name = "TestCreateCountTable"
    cds = CountDataSet(
        count_files=count_files,
        name=name,
        gene_column_name=gene_name,
        output_dir=tmpdir,
    )
    cds.create_count_table()
    expected_df = pd.DataFrame(
        {
            "barcode": {
                0: "AAAACCTCCCTGCCCAT",
                1: "ACGCAGACCCTCACTTT",
                2: "ACGCGCCAGACTTACGC",
                3: "AGAAGGAGAGCGAATAT",
                4: "AGATAACGAAACCACAC",
                5: "AGATCGCTGCTCGGGCG",
                6: "CTTCACTGTCATACGAA",
                7: "GAGCTAACCGATAACGG",
                8: "GATAGCTTGATGACGCA",
                9: "GTGACTCCGTCCAACAG",
                10: "TAAACAATGTACATAGA",
                11: "TATCGAACCACATCATA",
                12: "TCTGCAACGAGTTCAAC",
                13: "TGTTTTGGTAACGCTGC",
            },
            "Name": {
                0: "strB",
                1: "mobA",
                2: "traY",
                3: "hilC",
                4: "yjdB",
                5: "SL1344_4046",
                6: "garD",
                7: "SL1344_3745",
                8: "SL1344_3959",
                9: "hilD",
                10: "hilA",
                11: "hilD",
                12: "yhjV",
                13: "hilD",
            },
            "cf1": {
                0: 12571.0,
                1: 0.0,
                2: 19424.0,
                3: 0.0,
                4: 12611.0,
                5: 13661.0,
                6: 14052.0,
                7: 0.0,
                8: 14535.0,
                9: 0.0,
                10: 0.0,
                11: 18368.0,
                12: 15462.0,
                13: 26486.0,
            },
            "cf2": {
                0: 0.0,
                1: 13683.0,
                2: 18312.0,
                3: 0.0,
                4: 11722.0,
                5: 11964.0,
                6: 12879.0,
                7: 0.0,
                8: 12579.0,
                9: 0.0,
                10: 0.0,
                11: 17041.0,
                12: 15061.0,
                13: 23119.0,
            },
            "cf3": {
                0: 0.0,
                1: 0.0,
                2: 0.0,
                3: 21138.0,
                4: 0.0,
                5: 0.0,
                6: 20477.0,
                7: 17262.0,
                8: 18396.0,
                9: 28086.0,
                10: 19519.0,
                11: 51787.0,
                12: 22450.0,
                13: 59248.0,
            },
        }
    )
    cds.count_table.equals(expected_df)
