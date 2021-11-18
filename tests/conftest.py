import subprocess
import shlex
import pytest

@pytest.fixture
def parse_transposon():
    transposon = "GTGTATAAGAGACAG:17:13:before"
    tp2, bcLen, bc2tp2, before = transposon.split(":")
    return tp2,  int(bcLen), int(bc2tp2), before


@pytest.fixture
def get_structure(experiment='rbseq'):
    if experiment == 'rbseq':
        return 'GTGTATAAGAGACAG:17:13:before'
    elif experiment == 'wish':
        return 'GGAGGTTCACAATGTGGGAGGTCA:40:0:after'
    else:
        return None
#
# def capture(command_str):
#     command = shlex.split(command_str)
#     proc = subprocess.Popen(command, stdout=subprocess.PIPE,
#                             stderr=subprocess.PIPE,)
#     out, err = proc.communicate()
#     return out, err, proc
#
# @pytest.fixture
# def assert_files_are_same(file1, file2, verbose=False):
#     cmd_str = f'cmp {file1} {file2}'
#     out, err, proc = capture(cmd_str)
#     if verbose:
#         print(out)
#         print(err)
#     return proc.returncode