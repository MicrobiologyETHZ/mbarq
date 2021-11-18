from tnseq2.src.quantify import *
from tnseq2.src.sequence import stream_fa
from collections import Counter
import pandas as pd
import pickle
import subprocess
import shlex

# todo gives unclear error if provided with wrong transposon structure.
# todo when testing with wish data two tags don't get annotated properly ???

TESTDATA = "./tests/test_files"
EXPDATA = "./tests/expected_outcomes/counting"
OUTDIR = "./tests/tmp"


def capture(command_str):
    command = shlex.split(command_str)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,)
    out, err = proc.communicate()
    return out, err, proc


def assert_files_are_same(file1, file2, verbose=False):
    cmd_str = f'cmp {file1} {file2}'
    out, err, proc = capture(cmd_str)
    if verbose:
        print(out)
        print(err)
    assert proc.returncode == 0


def get_test_inserts(wish=True):
    if wish:
        return f'{TESTDATA}/LibraryA_sample_1.fq'
    return f'{TESTDATA}/dnaid2023_12_test.fasta'


def test_quantify_load_fq_barcodes():
    in_fq_file = get_test_inserts(wish=False)
    cnter = quantify_load_fq_barcodes(in_fq_file)
    with open(Path(EXPDATA)/'quntify_load_fq_barcodes.pkl', 'rb') as f:
        expected_cnter = Counter(pickle.load(f))
    assert cnter == expected_cnter



# def test_quantify_read_barcode_map_files():
#     outMap = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/testMap.tsv"
#
#     barcode_2_pos, barcode_2_abundance = quantify_barcodes.quantify_read_barcode_map_files(outMap)
#     with open('/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/B2P.json') as jf:
#         expectedB2P = json.load(jf)
#         expectedB2P = {key: tuple(val) for key, val in expectedB2P.items()}
#     with open('/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/B2A.json') as jf2:
#         expectedB2A = json.load(jf2)
#     assert barcode_2_pos == expectedB2P
#     assert barcode_2_abundance == expectedB2A
#
#
#
# def capture(command_str):
#     command = shlex.split(command_str)
#     proc = subprocess.Popen(command, stdout=subprocess.PIPE,
#                             stderr=subprocess.PIPE,)
#     out, err = proc.communicate()
#     return out, err, proc
#
# def to_str(bytes_or_str):
#     if isinstance(bytes_or_str, bytes):
#         value = bytes_or_str.decode('utf-8')
#     else:
#         value = bytes_or_str
#     return value
#
# def test_quantify_extract_annotated_correct():
#     with open('/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/B2P.json') as jf:
#         expectedB2P = json.load(jf)
#         barcode_2_position = {key: tuple(val) for key, val in expectedB2P.items()}
#     with open('/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/B2A.json') as jf2:
#         barcode_2_abundance = Counter(json.load(jf2))
#     print(barcode_2_abundance.most_common()[0:10])
#     in_fq_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/seqData/test100.fasta"
#     cnter = quantify_barcodes.quantify_load_fq_barcodes(in_fq_file)
#     max_edit_distance = 2
#     expectedOutFile = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/seqData/quantifyExtract.sorted.out"
#     testOutFile = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/quantifyExtractTest.out"
#     quantify_barcodes.quantify_extract_annotated_correct(barcode_2_position, barcode_2_abundance, cnter, testOutFile, max_edit_distance)
#     out, err, proc = capture(f'sort -k1 {testOutFile}')
#     with open(expectedOutFile, 'r') as o:
#         expectedLines = o.read()
#     assert to_str(out) == expectedLines
#
#
# def test_quantify():
#     in_fq_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/seqData/test100.fasta"
#     map_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/testMap.tsv"
#     testOutFile = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/quantifyQuantOut.tsv"
#     expectedOutFile = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/seqData/quantifyExtract.sorted.out"
#     tp2 = 'GTGTATAAGAGACAG'
#     bc2tp2 = 13
#     bcLen = 17
#     before = True
#     max_edit_distance = 2
#     quantify_barcodes.quantify(in_fq_file, map_file, testOutFile, tp2, bc2tp2, bcLen, before, max_edit_distance)
#     out, err, proc = capture(f'sort -k1 {testOutFile}')
#     with open(expectedOutFile, 'r') as o:
#         expectedLines = o.read()
#     assert to_str(out) == expectedLines
#
#


#with open(f'{EXPDATA}/quntify_load_fq_barcodes.pkl', 'wb') as fo:
#    pickle.dump(cnter, fo)