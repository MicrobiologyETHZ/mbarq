from mbarq.mapper import Mapper, AnnotatedMap
from collections import Counter
from pathlib import Path
import pandas as pd
import pickle

from test_utils import assert_files_are_same

EXPDATA = "./tests/expected_outcomes/mapping"
OUTDIR = "./tests/tmp"


def test_extract_barcodes_tn5(tn5_structure, map_test_data):
    """Check ...."""
    r1, genome = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome)
    seq_data.extract_barcodes()
    barcodes = [(bc.bc_seq, bc.host) for bc in seq_data.barcodes]
    with open(Path(EXPDATA)/'extract_barcode.pkl', 'rb') as f:
        expected_barcodes = pickle.load(f)
    assert Counter(barcodes) == Counter(expected_barcodes)


# todo test for _dereplicate_barcodes

def test__write_barcodes_to_fasta(tn5_structure, map_test_data, tmpdir):
    r1, genome = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir)
    seq_data.extract_barcodes()
    seq_data._write_barcodes_to_fasta()
    outFasta = tmpdir.join("library_13_1_1.fasta")
    expectedFasta = f'{EXPDATA}/prepare_extract_barcodes.fasta'
    assert_files_are_same(outFasta, expectedFasta)

# todo rename (expected) test files


def test__blast_barcode_host_genome(tn5_structure, map_test_data, tmpdir):
    global OUTDIR
    r1, genome = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir)
    seq_data.extract_barcodes()
    seq_data._write_barcodes_to_fasta()
    seq_data._blast_barcode_host()
    out_blastn_file = tmpdir.join("library_13_1_1.blastn")
    expected_blastn_file = f'{EXPDATA}/library_13_1_1.blastn'
    assert_files_are_same(out_blastn_file, expected_blastn_file)


def test__blast_barcode_host_genomedb(tn5_structure, map_test_data, tmpdir):
    r1, genome = map_test_data
    seq_data = Mapper(r1, tn5_structure, db=genome, output_dir=tmpdir)
    seq_data.extract_barcodes()
    seq_data._write_barcodes_to_fasta()
    seq_data._blast_barcode_host()
    out_blastn_file = tmpdir.join("library_13_1_1.blastn")
    expected_blastn_file = f'{EXPDATA}/library_13_1_1.blastn'
    assert_files_are_same(out_blastn_file, expected_blastn_file)


def test__find_most_likely_positions(tn5_structure, map_test_data, tmpdir):
    r1, genome = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=OUTDIR)
    seq_data.extract_barcodes()
    seq_data._write_barcodes_to_fasta()
    seq_data.temp_blastn_file = f'{EXPDATA}/library_13_1_1.blastn'
    seq_data._find_most_likely_positions(filter_below=0)
    expected_csv = f'{EXPDATA}/likely_positions.csv'
    tmp_csv = tmpdir.join("likely_positions.csv")
    seq_data.positions.to_csv(tmp_csv, index=False)
    assert_files_are_same(expected_csv, tmp_csv)


# todo add more test for edge cases (ex. multimappers)

def test_merge_colliding_bcs(tn5_structure, map_test_data, tmpdir):
    r1, genome = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir)
    seq_data.extract_barcodes()
    positions = pd.read_csv(f'{EXPDATA}/likely_positions.csv')
    seq_data.positions = positions
    seq_data._merge_colliding_barcodes()
    tmp_csv = tmpdir.join('merge_colliding_bcs.csv')
    seq_data.positions.to_csv(tmp_csv)
    expected_csv = f'{EXPDATA}/merge_colliding_bcs.csv'
    assert_files_are_same(tmp_csv, expected_csv)


def test_map(tn5_structure, map_test_data, tmpdir):
    r1, genome = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir)
    min_host_bases = 20
    filter_below = 0
    seq_data.map_insertions(min_host_bases, filter_below)
    out_map = tmpdir.join('library_13_1_1.map.csv')
    expected_map = f'{EXPDATA}/library_13_1_1.map.csv'
    assert_files_are_same(out_map, expected_map)


def test__find_annotation_overlaps_intersect(sl1344_gff, tmpdir):
    map_file = f'{EXPDATA}/library_13_1_1.map.csv'
    annotated_map = AnnotatedMap(map_file, sl1344_gff, "gene", ("ID", "Name", "locus_tag"),
                                 name='library_13_1_1', output_dir=tmpdir)
    annotated_map._find_annotation_overlaps(intersect=True)
    out_intersect = tmpdir.join('library_13_1_1.bed.intersect.tab')
    expected_intersect = f'{EXPDATA}/library_13_1_1.bed.intersect.tab'
    assert_files_are_same(out_intersect, expected_intersect)


def test__find_annotation_overlaps_closest(sl1344_gff, tmpdir):
    map_file = f'{EXPDATA}/library_13_1_1.map.csv'
    annotated_map = AnnotatedMap(map_file, sl1344_gff, "gene", ("ID", "Name", "locus_tag"),
                                 name='library_13_1_1', output_dir=tmpdir)
    annotated_map._find_annotation_overlaps(intersect=False)
    out_intersect = tmpdir.join('library_13_1_1.bed.intersect.tab')
    expected_intersect = f'{EXPDATA}/library_13_1_1.bed.closest.tab'
    assert_files_are_same(out_intersect, expected_intersect)


def test__add_bedtools_results_to_positions(sl1344_gff, tmpdir):
    map_file = f'{EXPDATA}/library_13_1_1.map.csv'
    annotated_map = AnnotatedMap(map_file, sl1344_gff, "gene", ("ID", "Name", "locus_tag"),
                                 name='library_13_1_1', output_dir=tmpdir)
    annotated_map._find_annotation_overlaps()
    annotated_map._add_bedtools_results_to_positions(intersect=True)
    out_map = tmpdir.join('library_13_1_1.annotated.csv')
    expected_map = f'{EXPDATA}/library_13_1_1.annotated.csv'
    assert_files_are_same(out_map, expected_map)


def test__add_bedtools_results_to_positions_with_closest(sl1344_gff, tmpdir):
    map_file = f'{EXPDATA}/library_13_1_1.map.csv'
    annotated_map = AnnotatedMap(map_file, sl1344_gff, "gene", ("ID", "Name", "locus_tag"),
                                 name='library_13_1_1_closest', output_dir=tmpdir)
    annotated_map._find_annotation_overlaps(intersect=False)
    annotated_map._add_bedtools_results_to_positions(intersect=False)
    out_map = tmpdir.join('library_13_1_1_closest.annotated.csv')
    expected_map = f'{EXPDATA}/library_13_1_1_closest.annotated.csv'
    assert_files_are_same(out_map, expected_map)

#
# def test_annotate_intersect(tmpdir):
#     r1, genome = get_test_file()
#     structure = get_structure()
#     gff_file = f'{TESTDATA}/ref/Salmonella_genome+plasmids.gff'
#     seq_data = Mapper(r1, structure, genome=genome, output_dir=tmpdir)
#     min_host_bases = 20
#     filter_below = 0
#     seq_data.map_insertions(min_host_bases, filter_below)
#     feature_type = 'gene'
#     identifiers = ('ID', 'Name', 'locus_tag')
#     annotated_map = AnnotatedMap(map_file=seq_data.map_file, annotation_file=gff_file, positions=seq_data.positions,
#                                  feature_type=feature_type, identifiers=identifiers, output_dir=tmpdir)
#     annotated_map.annotate(intersect=True)
#     out_map = tmpdir.join('library_13_1_1.map.annotated.csv')
#     expected_map = f'{EXPDATA}/library_13_1_1.annotated.csv'
#     assert_files_are_same(out_map, expected_map)
#
#
# # def test_annotate_intersect_no_gff(tmpdir):
# #     # todo convert this into a test
# #     gff_file = f'{TESTDATA}/ref/Salmonella_genome+plasmids.gffxx'
# #     feature_type = 'gene'
# #     identifiers = ('ID', 'Name', 'locus_tag')
# #     map_file = f'{EXPDATA}/library_13_1_1.map.csv'
# #     annotated_map = AnnotatedMap(map_file=map_file, annotation_file=gff_file,
# #                                  feature_type=feature_type, identifiers=identifiers, output_dir=tmpdir)
# #     annotated_map.annotate(intersect=True)
# #     #out_map = tmpdir.join('library_13_1_1.map.annotated.csv')
# #     #out_map = Path(OUTDIR)/'library_13_1_1.map.csv'
# #     #expected_map = f'{EXPDATA}/library_13_1_1.map.annotated.csv'
# #     #assert_files_are_same(out_map, expected_map)
# # test_annotate_intersect_no_gff(OUTDIR)
# #
# # #
# # def test_annotate_closest(tmpdir):
# #     r1, genome = get_test_file()
# #     structure = get_structure()
# #     gff_file = f'{TESTDATA}/ref/Salmonella_genome+plasmids.gff'
# #     seq_data = Mapper(r1, structure, genome=genome, output_dir=tmpdir)
# #     min_host_bases = 20
# #     filter_below = 0
# #     seq_data.map_insertions(min_host_bases, filter_below)
# #     feature_type = 'gene'
# #     identifiers = ('ID', 'Name', 'locus_tag')
# #     annotated_map = AnnotatedMap(map_file=seq_data.map_file, annotation_file=gff_file, positions=seq_data.positions,
# #                                  feature_type=feature_type, identifiers=identifiers, output_dir=tmpdir)
# #     annotated_map.annotate(intersect=False)
# #     out_map = tmpdir.join('library_13_1_1.map.annotated.csv')
# #     expected_map = f'{EXPDATA}/library_13_1_1.map.closest.annotated.csv'
# #     assert_files_are_same(out_map, expected_map)
