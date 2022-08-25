import subprocess
import shlex
import pytest
from pathlib import Path

#root = Path("/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/Projects_NCCR/ref/mbarq_test_data")
root = Path("/Users/ansintsova/Downloads/test_mbarq/mbarq_test_data")
@pytest.fixture
def tn5_structure():
    return 'B17N13GTGTATAAGAGACAG'

@pytest.fixture
def wish_structure():
    return 'GGAGGTTCACAATGTGGGAGGTCAB40'


@pytest.fixture
def map_test_data():
    seq_file = root/"dnaid1315/test_data/library_11_1_FKDL202598974-1a-D701-AK1682_HHG5YDSXY_L4_subsample_1.fq.gz"
    genome_file = root/"dnaid1315/ref/GCA_000210855.2_ASM21085v2_genomic.fna"
    return (seq_file, genome_file)


@pytest.fixture
def sl1344_gff():
    return root/"dnaid1315/ref/GCA_000210855.2_ASM21085v2_genomic.gff"


@pytest.fixture
def count_test_data_tn5():
    small_count_file = root/"dnaid1315/test_data/count_test.fasta.gz"
    seq_file = root/"dnaid1315/test_data/dnaid1315_124_subsample.fasta.gz"
    map_file = root/"dnaid1315/ref/library_11_1.annotated.csv"
    return (small_count_file, seq_file, map_file)


@pytest.fixture
def count_test_data_wish():
    seq_file = f'./tests/test_files/LibraryA_pilot2.fq.gz'
    map_file = f'./tests/test_files/20210520_BarcodeList.csv'
    return (seq_file, map_file)


@pytest.fixture
def merge_test_data_tn5():
    count_files = [f for f in (root/"dnaid1315/ref").glob("dnaid1315_*_mbarq_counts.csv")]
    return count_files, 'Name'

@pytest.fixture
def analysis_test_data_tn5():
    control1col = root/"dnaid1315/ref/controls_1col.csv"
    control2col = root / "dnaid1315/ref/controls_2col.csv"
    control3col = root / "dnaid1315/ref/controls_3col.csv"
    control2col_short = root / "dnaid1315/ref/controls_2col_short.csv"
    merged_counts = root/"dnaid1315/ref/example_mbarq_merged_counts.csv"
    sample_data = root/"dnaid1315/ref/example_sample_data.csv"
    no_wt = root/"dnaid1315/ref/controls_3col_no_wt.csv"
    return control1col,control2col,control3col,control2col_short,merged_counts,sample_data, no_wt

@pytest.fixture
def dnaid1315_expected_outcomes():
    return root/"dnaid1315/expected_outcomes"

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