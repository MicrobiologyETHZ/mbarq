import subprocess
import shlex
from pathlib import Path
from tnseq2.src.sequence import stream_fa


TESTDATA = "./tests/test_files"
EXPDATA = "./tests/expected_outcomes"
OUTDIR = "./tests/tmp"

# todo not DRY, need to refactor


def capture(command_str):
    command = shlex.split(command_str)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,)
    out, err = proc.communicate()
    return out, err, proc


def to_str(bytes_or_str):
    if isinstance(bytes_or_str, bytes):
        value = bytes_or_str.decode('utf-8')
    else:
        value = bytes_or_str
    return value


def assert_files_are_same(file1, file2, verbose=False):
    cmd_str = f'cmp {file1} {file2}'
    out, err, proc = capture(cmd_str)
    if verbose:
        print(out)
        print(err)
    assert proc.returncode == 0


def get_test_inserts(files=False):
    r1 = f'{TESTDATA}/library_13_1_1.fq'
    r2 = f'{TESTDATA}/library_13_1_2.fq'
    if files:
        return r1, r2
    inserts = zip(stream_fa(r1), stream_fa(r2))
    return inserts


def get_ref_data_map():
    genome = f'{TESTDATA}/ref/Salmonella_genome_FQ312003.1_SL1344.fasta'
    gff_file = f'{TESTDATA}/ref/Salmonella_genome+plasmids.gff'
    return genome, gff_file


# def get_test_data_quant():
#     r = "/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/DEV/TNSEQ_DEV/data/processed/tnseq/tnseq_pipeline/input_quantify/1315-107-library11_1-TV3371B-inoculum.fasta.gz"
#     outMap = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/testMap.tsv"
#     expectedQuant = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/outQuant.tsv"
#     outQuant = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/testQuant.tsv"
#     return r, outMap, expectedQuant, outQuant


def test_cli_map():
    r1, r2 = get_test_inserts(files=True)
    genome, gff_file = get_ref_data_map()
    cmd_str = f'tnseq2 maplib -f {r1} -r {r2} -a {gff_file} ' \
              f' -o {OUTDIR} -g {genome} -n Test1 -l 0'
    _, err, proc = capture(cmd_str)
    if proc.returncode != 0:
        print(to_str(err))
    assert proc.returncode == 0
    outMap = Path(OUTDIR)/"Test1.barcode_map.annotated.csv"
    expectedMap = Path(EXPDATA)/"mapping/merge_colliding_bcs.annotated.csv"
    assert_files_are_same(outMap, expectedMap)




#
# def test_cli_quantify():
#
#     r, outMap, expectedQuant, outQuant = get_test_data_quant()
#     cmd_str = f'tnseq quantify -r {r}  -o {outQuant} -b {outMap} -f'
#     _, err, proc = capture(cmd_str)
#     assert proc.returncode == 0
#     out, err, proc = capture(f'sort -k1 {outQuant}')
#     with open(expectedQuant, 'r') as o:
#        expectedLines = o.read()
#     assert to_str(out) == expectedLines
