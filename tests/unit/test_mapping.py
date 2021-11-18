from tnseq2.src.mapping import *
from tnseq2.src.sequence import stream_fa
from collections import Counter
import pandas as pd
import pickle
import subprocess
import shlex

TESTDATA = "./tests/test_files"
EXPDATA = "./tests/expected_outcomes/mapping"
OUTDIR = "./tests/tmp"

# todo add pytest tmp directory, will help with teardown
# todo generalize inputs between functions

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


def get_test_inserts(files=False):
    r1 = f'{TESTDATA}/library_13_1_1.fq'
    r2 = f'{TESTDATA}/library_13_1_2.fq'
    if files:
        return r1, r2
    inserts = zip(stream_fa(r1), stream_fa(r2))
    return inserts


def parse_transposon(transposon="GTGTATAAGAGACAG:17:13:before"):
    tp2, bcLen, bc2tp2, before = transposon.split(":")
    return tp2,  int(bcLen), int(bc2tp2), before


def test_extract_barcodes():
    """Check ...."""
    inserts = get_test_inserts()
    tp2,  bcLen, bc2tp2, before = parse_transposon()
    barcodes = extract_barcodes(inserts, tp2, bc2tp2, bcLen, before)
    with open(Path(EXPDATA)/'extract_barcode.pkl', 'rb') as f:
        expected_barcodes = pickle.load(f)
    assert Counter(barcodes) == Counter(expected_barcodes)


def test_prepare_extract_barcodes(tmpdir):
    r1, r2 = get_test_inserts(files=True)
    outFasta = tmpdir.join('prepare_extract_barcodes.fasta')
    prepare_extract_barcodes(r1, r2, outFasta)
    expectedFasta = f'{EXPDATA}/prepare_extract_barcodes.fasta'
    assert_files_are_same(outFasta, expectedFasta)


def test_map_host_map(tmpdir):
    temp_fasta_file = f'{EXPDATA}/prepare_extract_barcodes.fasta'
    temp_blastn_file = tmpdir.join('map_host_map.blast')
    genome = f'{TESTDATA}/ref/Salmonella_genome_FQ312003.1_SL1344.fasta'
    blastdb = ''
    db_rc, blast_rc = map_host_map(temp_fasta_file, temp_blastn_file, genome, blastdb)
    assert db_rc == 0
    assert blast_rc == 0
    blastdb = f'{TESTDATA}/ref/Salmonella_genome_FQ312003.1_SL1344.fasta'
    db_rc, blast_rc = map_host_map(temp_fasta_file, temp_blastn_file, genome, blastdb)
    assert db_rc == 0
    assert blast_rc == 0
    expected_blastn_file = f'{EXPDATA}/map_host_map.blast'
    assert_files_are_same(temp_blastn_file, expected_blastn_file)


def test_new_map_annotate(tmpdir):
    blast_file = f'{EXPDATA}/map_host_map.blast'
    tmp_df = new_map_annotate(blast_file, filter_below=0)
    expected_csv = f'{EXPDATA}/map_annotate.csv'
    tmp_csv = tmpdir.join("map_annotate.csv")
    tmp_df.to_csv(tmp_csv)
    assert_files_are_same(expected_csv, tmp_csv)

# todo add more test for edge cases (ex. multimappers)


def test_merge_colliding_bcs(tmpdir):
    pot_pos = pd.read_csv(f'{EXPDATA}/map_annotate.csv', index_col=0)
    tmp_df = merge_colliding_bcs(pot_pos, edit_dist=3)
    tmp_csv = tmpdir.join('merge_colliding_bcs.csv')
    tmp_df.to_csv(tmp_csv)
    expected_csv = f'{EXPDATA}/merge_colliding_bcs.csv'
    assert_files_are_same(tmp_csv, expected_csv)


def test_add_gff_annotations(tmpdir):
    barcode_map = pd.read_csv(f'{EXPDATA}/merge_colliding_bcs.csv', index_col=0)
    bed_file = tmpdir.join('tmp_file.bed')
    gff_file = f'{TESTDATA}/ref/Salmonella_genome+plasmids.gff'
    tmp_output_map = tmpdir.join('add_gff_annotations.tab')
    return_code = add_gff_annotations(barcode_map, bed_file, gff_file, tmp_output_map)
    assert return_code == 0
    expected_output_map = f'{EXPDATA}/add_gff_annotations.tab'
    assert_files_are_same(tmp_output_map, expected_output_map)


def test_add_gene_annotations(tmpdir):
    bedfile = f'{EXPDATA}/add_gff_annotations.tab'
    barcode_map = pd.read_csv(f'{EXPDATA}/merge_colliding_bcs.csv', index_col=0)
    tmp_barcode_file = tmpdir.join('merge_colliding_bcs.csv')
    add_gene_annotations(bedfile, barcode_map, tmp_barcode_file)
    tmp_final_file = tmpdir.join('merge_colliding_bcs.annotated.csv')
    expected_final_file = f'{EXPDATA}/merge_colliding_bcs.annotated.csv'
    assert_files_are_same(tmp_final_file, expected_final_file)


def test_map_host_new(tmpdir):
    fasta_file = f'{EXPDATA}/prepare_extract_barcodes.fasta'
    blastn_file = tmpdir.join('map_host_map.blast')
    genome = f'{TESTDATA}/ref/Salmonella_genome_FQ312003.1_SL1344.fasta'
    blast_threads = 1
    tmp_bed_file = tmpdir.join('tmp_file.bed')
    output_map = tmpdir.join('add_gff_annotations.tab')
    barcode_file = tmpdir.join('merge_colliding_bcs.csv')
    gff_file = f'{TESTDATA}/ref/Salmonella_genome+plasmids.gff'
    filter_below = 0
    map_host_new(fasta_file, blastn_file, genome,
                 blast_threads, tmp_bed_file,
                 output_map, barcode_file, gff_file, filter_below)
    expected_barcode_file = f'{EXPDATA}/merge_colliding_bcs.csv'
    expected_annotated_barcode_file = f'{EXPDATA}/merge_colliding_bcs.annotated.csv'
    assert_files_are_same(barcode_file, expected_barcode_file)
    assert_files_are_same(Path(barcode_file).with_suffix('.annotated.csv'), expected_annotated_barcode_file)


def test_map(tmpdir):
    r1_file, r2_file = get_test_inserts(files=True)
    name='Test'
    output_dir = tmpdir.mkdir('tmp')
    transposon = "GTGTATAAGAGACAG:17:13:before"
    genome = f'{TESTDATA}/ref/Salmonella_genome_FQ312003.1_SL1344.fasta'
    gff_file = f'{TESTDATA}/ref/Salmonella_genome+plasmids.gff'
    map(r1_file, r2_file, name, output_dir, transposon, genome, gff_file,
        min_host_bases=20, blast_threads=1, filter_below=0)
    expected_barcode_file = f'{EXPDATA}/merge_colliding_bcs.csv'
    expected_annotated_barcode_file = f'{EXPDATA}/merge_colliding_bcs.annotated.csv'
    assert_files_are_same(Path(output_dir)/'Test.barcode_map.csv', expected_barcode_file)
    assert_files_are_same(Path(output_dir)/'Test.barcode_map.annotated.csv', expected_annotated_barcode_file)



