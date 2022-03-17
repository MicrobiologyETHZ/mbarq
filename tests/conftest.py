import subprocess
import shlex
import pytest

@pytest.fixture
def tn5_structure():
    return 'B17N13GTGTATAAGAGACAG'

@pytest.fixture
def wish_structure():
    return 'GGAGGTTCACAATGTGGGAGGTCAB40'


@pytest.fixture
def map_test_data():
    seq_file = "./tests/test_files/library_13_1_1.fq.gz"
    genome_file = "./tests/test_files/GCA_000210855.2_ASM21085v2_genomic.fna"
    return (seq_file, genome_file)


@pytest.fixture
def sl1344_gff():
    return "./tests/test_files/GCA_000210855.2_ASM21085v2_genomic.gff"


@pytest.fixture
def count_test_data_tn5():
    seq_file = "./tests/test_files/dnaid2023_12_test.fasta"
    map_file = "/Users/ansintsova/git_repos/mbarq/tests/test_files/library_13_1.barcode_map.annotated.csv"
    return (seq_file, map_file)


@pytest.fixture
def count_test_data_wish():
    seq_file = f'./tests/test_files/LibraryA_pilot2.fq.gz'
    map_file = f'./tests/test_files/20210520_BarcodeList.csv'
    return (seq_file, map_file)


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