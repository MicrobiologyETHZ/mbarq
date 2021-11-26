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
def test__find_annotation_overlaps(tmpdir):

    r1, genome = get_test_file()
    structure = get_structure()
    gff_file = f'{TESTDATA}/ref/Salmonella_genome+plasmids.gff'
    seq_data = Mapper(r1, structure, genome=genome, annotation_file=gff_file, output_dir=tmpdir)
    seq_data.positions = pd.read_csv(f'{EXPDATA}/merge_colliding_bcs_new.csv', index_col=0)
    seq_data._find_annotation_overlaps()
    out_tab = tmpdir.join('library_13_1_1.bed.intersect.tab')
    expected_tab = f'{EXPDATA}/library_13_1_1.bed.intersect.tab'
    assert_files_are_same(out_tab, expected_tab)

def test__find_closest_feature(tmpdir):
    r1, genome = get_test_file()
    structure = get_structure()
    gff_file = f'{TESTDATA}/ref/Salmonella_genome+plasmids.gff'
    seq_data = Mapper(r1, structure, genome=genome, annotation_file=gff_file, output_dir=tmpdir)
    seq_data.positions = pd.read_csv(f'{EXPDATA}/merge_colliding_bcs_new.csv', index_col=0)
    seq_data._find_closest_feature()
    out_tab = tmpdir.join('library_13_1_1.bed.closest.tab')
    expected_tab = f'{EXPDATA}/library_13_1_1.bed.closest.tab'
    assert_files_are_same(out_tab, expected_tab)


def test__add_bedintersect_results_to_positions(tmpdir):
    global OUTDIR
    r1, genome = get_test_file()
    structure = get_structure()
    seq_data = Mapper(r1, structure, genome=genome, output_dir=tmpdir)
    seq_data.identifiers = ('ID', 'Name', 'locus_tag')
    seq_data.temp_bed_results_file = f'{EXPDATA}/library_13_1_1.bed.intersect.tab'
    seq_data.positions = pd.read_csv(f'{EXPDATA}/merge_colliding_bcs_new.csv', index_col=0)
    seq_data._add_bedintersect_results_to_positions('gene')
    seq_data.write_mapfile()
    out_map = tmpdir.join('library_13_1_1.map.annotated.csv')
    expected_map = f'{EXPDATA}/library_13_1_1.map.annotated.csv'
    assert_files_are_same(out_map, expected_map)


def test__add_bedintersect_results_to_positions_with_closest(tmpdir):

    r1, genome = get_test_file()
    structure = get_structure()
    seq_data = Mapper(r1, structure, genome=genome, output_dir=tmpdir)
    seq_data.identifiers = ('ID', 'Name', 'locus_tag')
    seq_data.temp_bed_closest_file = f'{EXPDATA}/library_13_1_1.bed.closest.tab'
    seq_data.positions = pd.read_csv(f'{EXPDATA}/merge_colliding_bcs_new.csv', index_col=0)
    seq_data._add_bedintersect_results_to_positions('gene', False)
    seq_data.write_mapfile()
    out_map = tmpdir.join('library_13_1_1.map.annotated.csv')
    expected_map = f'{EXPDATA}/library_13_1_1.map.closest.annotated.csv'
    assert_files_are_same(out_map, expected_map)


def test_map(tmpdir):
    global OUTDIR
    r1, genome = get_test_file()
    structure = get_structure()
    gff_file = f'{TESTDATA}/ref/Salmonella_genome+plasmids.gff'
    seq_data = Mapper(r1, structure, genome=genome, annotation_file=gff_file, output_dir=tmpdir)
    min_host_bases = 20
    filter_below = 0
    seq_data.map_insertions(min_host_bases, filter_below)
    seq_data.write_mapfile()
    out_map = tmpdir.join('library_13_1_1.map.csv')
    expected_map = f'{EXPDATA}/library_13_1_1.map.csv'
    assert_files_are_same(out_map, expected_map)


def test_annotate_intersect(tmpdir):
    r1, genome = get_test_file()
    structure = get_structure()
    gff_file = f'{TESTDATA}/ref/Salmonella_genome+plasmids.gff'
    seq_data = Mapper(r1, structure, genome=genome, output_dir=tmpdir)
    min_host_bases = 20
    filter_below = 0
    seq_data.map_insertions(min_host_bases, filter_below)
    feature_type = 'gene'
    identifiers = ('ID', 'Name', 'locus_tag')
    seq_data.annotate(annotations_file=gff_file, feature_type=feature_type,
                      identifiers=identifiers, intersect=True)
    seq_data.write_mapfile()
    out_map = tmpdir.join('library_13_1_1.map.annotated.csv')
    expected_map = f'{EXPDATA}/library_13_1_1.map.annotated.csv'
    assert_files_are_same(out_map, expected_map)

#
def test_annotate_closest(tmpdir):
    r1, genome = get_test_file()
    structure = get_structure()
    gff_file = f'{TESTDATA}/ref/Salmonella_genome+plasmids.gff'
    seq_data = Mapper(r1, structure, genome=genome, output_dir=tmpdir)
    min_host_bases = 20
    filter_below = 0
    seq_data.map_insertions(min_host_bases, filter_below)
    feature_type = 'gene'
    identifiers = ('ID', 'Name', 'locus_tag')
    seq_data.annotate(annotations_file=gff_file, feature_type=feature_type,
                      identifiers=identifiers, intersect=False)
    seq_data.write_mapfile()
    out_map = tmpdir.join('library_13_1_1.map.annotated.csv')
    expected_map = f'{EXPDATA}/library_13_1_1.map.closest.annotated.csv'
    assert_files_are_same(out_map, expected_map)
