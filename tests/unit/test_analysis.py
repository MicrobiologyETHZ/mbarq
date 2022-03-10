from mbarq.analysis import CountDataSet
from collections import Counter
from pathlib import Path
import pandas as pd
import pickle
import subprocess
import shlex
import json
import random



TMPDIR = "/Users/ansintsova/git_repos/mbarq/tests/tmp"


def get_test_file(experiment='rbseq'):
    if experiment == 'rbseq':
        return [f'{TMPDIR}/dnaid2023_12_test2_mbarq_counts.csv', f'{TMPDIR}/dnaid2023_12_test_mbarq_counts.csv']
    else:
        return None

def test__merge_count_files():
    count_files = get_test_file()
    name = 'TestingMerge'
    gene_name = 'ShortName'
    cds = CountDataSet(count_files=count_files, name=name, gene_column_name=gene_name, output_dir=TMPDIR)
    fdf = cds._merge_count_files()
    print(fdf.head())


test__merge_count_files()