import subprocess
import shlex
import pytest
from pathlib import Path
from mbarq.core import Barcode

root = Path("tests/mbarq_test_data")

@pytest.fixture
def tn5_structure():
    return 'B17N13GTGTATAAGAGACAG'

@pytest.fixture
def wish_structure():
    return 'GGAGGTTCACAATGTGGGAGGTCAB40'


@pytest.fixture
def map_test_data():
    seq_file = root/"mapping/mapping_reads.fastq.gz"
    genome_file = root/"mapping/test_genome.fna"
    barcodes = [Barcode('B17N13GTGTATAAGAGACAG', sequence, host) for sequence, host in [('ATGAAAAGTATGAATCC',
                                                                                        'GAATAGACATTTTAACACTCCGAAATCATTATAAATGAATAATTAACACAAGGAGTGTAAGCGCCCTGTCAGGAGGGTAAAAATCAAAGCACATCATTTAAAACC'),
                                                                                       ('ATGAAAAGTATGAATCC',
                                                                                        'GAATAGACATTTTAACACTCCGAAATCATTATAAATGAATAATTAACACAAGGAGTGTAAGCGCCCTGTCAGGAGGGTAAAAATCAAAGCACATCATTTAAAACC'),
                                                                                       ('ATGAAAAGTATGAATCC',
                                                                                        'GAATAGACATTTTAACACTCCGAAATCATTATAAATGAATAATTAACACAAGGAGTGTAAGCGCCCTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATTA'),
                                                                                       ('TTGGGATCCCACCATTT',
                                                                                        'GAGTAAGCGGGCGCTGAGAGGTGTTGTTTTCTCTTCGTTAGACGGTGTTGTTAACCTCATTTTTATGATTTTTATATCATCTAAAAAGATGATGTTTTGTGATTA'),
                                                                                       ('CGTCAGGGCAGCGAACA',
                                                                                        'AACTCAGTCTAACGCCAAGGGTCTGCTGGGCGCGCTGCGTGATATGCAGGCAAAAGCGAAAGCCGCAGGTCACACGCTGGCGCTCTCCATGAGTATCGGCGGCTG'),
                                                                                       ('CGTCAGGGCAGCGAACA',
                                                                                        'AACTCAGTCTAACGCCAAGGGTCTGCTGGGCGCGCTGCGTGATATGCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATTACTCGATCTCGTATGCCGTCT'),
                                                                                       ]]
    mapping_dir = root/'mapping'
    return seq_file, genome_file, barcodes, mapping_dir


@pytest.fixture
def test_gff():
    return root/"mapping/map_annotate.gff"


@pytest.fixture
def count_test_data():
    seq_file = root/"counting/dnaid1315_124_subsample.fasta.gz"
    map_file = root/"counting/library.annotated.csv"
    return (seq_file, map_file)


# @pytest.fixture
# def count_test_data_wish():
#     seq_file = f'./tests/test_files/LibraryA_pilot2.fq.gz'
#     map_file = f'./tests/test_files/20210520_BarcodeList.csv'
#     return (seq_file, map_file)


@pytest.fixture
def merge_test_data():
    count_files = [f for f in (root/"counting").glob("*_mbarq_counts.csv")]
    return count_files, 'Name'

@pytest.fixture
def analysis_test_data():
    control1col = root/"analysis/controls_1col.csv.gz"
    control2col = root / "analysis/controls_2col.csv.gz"
    control3col = root / "analysis/controls_3col.csv.gz"
    control2col_short = root / "analysis/controls_2col_short.csv.gz"
    merged_counts = root/"analysis/example_mbarq_merged_counts.csv.gz"
    sample_data = root/"analysis/example_sample_data.csv"
    no_wt = root/"analysis/controls_3col_no_wt.csv.gz"
    control3bc = root/"analysis/controls_3bc.csv.gz"
    return control1col,control2col,control3col,control2col_short,merged_counts,sample_data, no_wt, control3bc


def capture(command_str):
    command = shlex.split(command_str)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,)
    out, err = proc.communicate()
    return out, err, proc
#
# @pytest.fixture
# def assert_files_are_same(file1, file2, verbose=False):
#     cmd_str = f'cmp {file1} {file2}'
#     out, err, proc = capture(cmd_str)
#     if verbose:
#         print(out)
#         print(err)
#     return proc.returncode