from mbarq.mapper import Mapper
from collections import Counter
from pathlib import Path
import pandas as pd
import pickle
import subprocess
import shlex

TESTDATA = "./tests/test_files"
EXPDATA = "./tests/expected_outcomes/mapping"
OUTDIR = "./tests/tmp"

# todo add pytest tmp directory, will help with teardown
# todo generalize inputs between functions
# todo add common functions to conftest.py


def get_structure(experiment='rbseq'):
    if experiment == 'rbseq':
        return 'GTGTATAAGAGACAG:17:13:before'
    elif experiment == 'wish':
        return 'GGAGGTTCACAATGTGGGAGGTCA:40:0:after'
    else:
        return None


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


def get_test_file(experiment = 'rbseq'):
    if experiment == 'rbseq':
        return f'{TESTDATA}/library_13_1_1.fq', f'{TESTDATA}/ref/Salmonella_genome_FQ312003.1_SL1344.fasta'
    else:
        return None

#
# def get_test_inserts(files=False, experiment='rbseq'):
#     r1 = f'{TESTDATA}/library_13_1_1.fq'
#     structure = get_structure(experiment)
#     if files:
#         return r1
#     seq_data = Mapper(r1, structure)
#     inserts = seq_data.stream_seq_file()
#     return inserts


def test_extract_barcodes_rbseq():
    """Check ...."""
    r1, genome = get_test_file()
    structure = get_structure()
    seq_data = Mapper(r1, structure, genome=genome)
    seq_data.extract_barcodes()
    barcodes = [(bc.bc_seq, bc.host) for bc in seq_data.barcodes]
    with open(Path(EXPDATA)/'extract_barcode.pkl', 'rb') as f:
        expected_barcodes = pickle.load(f)
    assert Counter(barcodes) == Counter(expected_barcodes)


# todo test for _dereplicate_barcodes

def test__write_barcodes_to_fasta(tmpdir):
    r1, genome = get_test_file()
    structure = get_structure()
    seq_data = Mapper(r1, structure, genome=genome, output_dir=tmpdir)
    seq_data.extract_barcodes()
    seq_data._write_barcodes_to_fasta()
    outFasta = tmpdir.join("library_13_1_1.fasta")
    expectedFasta = f'{EXPDATA}/prepare_extract_barcodes.fasta'
    assert_files_are_same(outFasta, expectedFasta)

# todo rename (expected) test files


def test__blast_barcode_host_genome(tmpdir):
    global OUTDIR
    r1, genome = get_test_file()
    structure = get_structure()
    seq_data = Mapper(r1, structure, genome=genome, output_dir=tmpdir)
    seq_data.extract_barcodes()
    seq_data._write_barcodes_to_fasta()
    seq_data._blast_barcode_host()
    out_blastn_file = tmpdir.join("library_13_1_1.blastn")
    expected_blastn_file = f'{EXPDATA}/map_host_map.blast'
    assert_files_are_same(out_blastn_file, expected_blastn_file)


def test__blast_barcode_host_genomedb(tmpdir):
    r1, genome = get_test_file()
    structure = get_structure()
    blastdb = f'{TESTDATA}/ref/Salmonella_genome_FQ312003.1_SL1344.fasta'
    seq_data = Mapper(r1, structure, db=blastdb, output_dir=tmpdir)
    seq_data.extract_barcodes()
    seq_data._write_barcodes_to_fasta()
    seq_data._blast_barcode_host()
    out_blastn_file = tmpdir.join("library_13_1_1.blastn")
    expected_blastn_file = f'{EXPDATA}/map_host_map.blast'
    assert_files_are_same(out_blastn_file, expected_blastn_file)


def test__find_most_likely_positions(tmpdir):
    #blast_file = f'{EXPDATA}/map_host_map.blast'

    r1, genome = get_test_file()
    structure = get_structure()
    blastdb = f'{TESTDATA}/ref/Salmonella_genome_FQ312003.1_SL1344.fasta'
    seq_data = Mapper(r1, structure, db=blastdb, output_dir=tmpdir)
    seq_data.extract_barcodes()
    seq_data._write_barcodes_to_fasta()
    seq_data._blast_barcode_host()
    seq_data._find_most_likely_positions(0)
    #tmp_df = new_map_annotate(blast_file, filter_below=0)
    expected_csv = f'{EXPDATA}/map_annotate.csv'
    tmp_csv = tmpdir.join("map_annotate.csv")
    seq_data.positions.to_csv(tmp_csv)
    assert_files_are_same(expected_csv, tmp_csv)

# todo add more test for edge cases (ex. multimappers)
# todo will have to re-write this after change filtering and insertion sites


def test_merge_colliding_bcs(tmpdir):
    r1, genome = get_test_file()
    structure = get_structure()
    seq_data = Mapper(r1, structure, genome=genome, output_dir=tmpdir)
    seq_data.extract_barcodes()
    seq_data._write_barcodes_to_fasta()
    seq_data._blast_barcode_host()
    seq_data._find_most_likely_positions(0)
    seq_data._merge_colliding_barcodes()
    tmp_csv = tmpdir.join('merge_colliding_bcs.csv')
    seq_data.positions.to_csv(tmp_csv)
    expected_csv = f'{EXPDATA}/merge_colliding_bcs_new.csv'
    assert_files_are_same(tmp_csv, expected_csv)
#
#
# def test_add_gff_annotations(tmpdir):
#     barcode_map = pd.read_csv(f'{EXPDATA}/merge_colliding_bcs.csv', index_col=0)
#     bed_file = tmpdir.join('tmp_file.bed')
#     gff_file = f'{TESTDATA}/ref/Salmonella_genome+plasmids.gff'
#     tmp_output_map = tmpdir.join('add_gff_annotations.tab')
#     return_code = add_gff_annotations(barcode_map, bed_file, gff_file, tmp_output_map)
#     assert return_code == 0
#     expected_output_map = f'{EXPDATA}/add_gff_annotations.tab'
#     assert_files_are_same(tmp_output_map, expected_output_map)
#
#
# def test_add_gene_annotations(tmpdir):
#     bedfile = f'{EXPDATA}/add_gff_annotations.tab'
#     barcode_map = pd.read_csv(f'{EXPDATA}/merge_colliding_bcs.csv', index_col=0)
#     tmp_barcode_file = tmpdir.join('merge_colliding_bcs.csv')
#     add_gene_annotations(bedfile, barcode_map, tmp_barcode_file)
#     tmp_final_file = tmpdir.join('merge_colliding_bcs.annotated.csv')
#     expected_final_file = f'{EXPDATA}/merge_colliding_bcs.annotated.csv'
#     assert_files_are_same(tmp_final_file, expected_final_file)
#
#

#
# def test_map(tmpdir):
#     r1_file, r2_file = get_test_inserts(files=True)
#     name='Test'
#     output_dir = tmpdir.mkdir('tmp')
#     transposon = "GTGTATAAGAGACAG:17:13:before"
#     genome = f'{TESTDATA}/ref/Salmonella_genome_FQ312003.1_SL1344.fasta'
#     gff_file = f'{TESTDATA}/ref/Salmonella_genome+plasmids.gff'
#     map(r1_file, r2_file, name, output_dir, transposon, genome, gff_file,
#         min_host_bases=20, blast_threads=1, filter_below=0)
#     expected_barcode_file = f'{EXPDATA}/merge_colliding_bcs.csv'
#     expected_annotated_barcode_file = f'{EXPDATA}/merge_colliding_bcs.annotated.csv'
#     assert_files_are_same(Path(output_dir)/'Test.barcode_map.csv', expected_barcode_file)
#     assert_files_are_same(Path(output_dir)/'Test.barcode_map.annotated.csv', expected_annotated_barcode_file)



