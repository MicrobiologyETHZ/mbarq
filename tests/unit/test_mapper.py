from mbarq.mapper import Mapper, AnnotatedMap
from collections import Counter
from pathlib import Path
import pandas as pd
import pickle
import pretty_errors
from test_utils import assert_files_are_same


def test_extract_barcodes_tn5(tn5_structure, map_test_data, dnaid1315_expected_outcomes, tmpdir):
    """Check ...."""
    r1, genome = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir)
    seq_data.extract_barcodes()
    barcodes = [(bc.bc_seq, bc.host) for bc in seq_data.barcodes]
    expected_file = dnaid1315_expected_outcomes/"extract_barcodes_tn5.pkl"
    with open(expected_file, 'rb') as f:
        expected_barcodes = pickle.load(f)
    assert Counter(barcodes) == Counter(expected_barcodes)


# todo test for _dereplicate_barcodes

def test__write_barcodes_to_fasta(tn5_structure, map_test_data, tmpdir, dnaid1315_expected_outcomes):
    r1, genome = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir)
    seq_data.extract_barcodes()
    seq_data._write_barcodes_to_fasta()
    outFasta = tmpdir.join("library_11_1_FKDL202598974-1a-D701-AK1682_HHG5YDSXY_L4_subsample_1.fasta")
    expectedFasta = f'{dnaid1315_expected_outcomes}/library_11_1_FKDL202598974-1a-D701-AK1682_HHG5YDSXY_L4_subsample_1.fasta'
    assert_files_are_same(outFasta, expectedFasta)


def test__blast_barcode_host_genome(tn5_structure, map_test_data, tmpdir, dnaid1315_expected_outcomes):
    r1, genome = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir)
    seq_data.extract_barcodes()
    seq_data._write_barcodes_to_fasta()
    seq_data._blast_barcode_host()
    out_blastn_file = tmpdir.join("library_11_1_FKDL202598974-1a-D701-AK1682_HHG5YDSXY_L4_subsample_1.blastn")
    expected_blastn_file = dnaid1315_expected_outcomes/"library_11_1_FKDL202598974-1a-D701-AK1682_HHG5YDSXY_L4_subsample_1.blastn"
    assert_files_are_same(out_blastn_file, expected_blastn_file)


def test__blast_barcode_host_genomedb(tn5_structure, map_test_data, tmpdir, dnaid1315_expected_outcomes):
    r1, genome = map_test_data
    seq_data = Mapper(r1, tn5_structure, db=genome, output_dir=tmpdir)
    seq_data.extract_barcodes()
    seq_data._write_barcodes_to_fasta()
    seq_data._blast_barcode_host()
    out_blastn_file = tmpdir.join("library_11_1_FKDL202598974-1a-D701-AK1682_HHG5YDSXY_L4_subsample_1.blastn")
    expected_blastn_file = dnaid1315_expected_outcomes/"library_11_1_FKDL202598974-1a-D701-AK1682_HHG5YDSXY_L4_subsample_1.blastn"
    assert_files_are_same(out_blastn_file, expected_blastn_file)


def test__find_most_likely_positions(tn5_structure, map_test_data, tmpdir, dnaid1315_expected_outcomes):
    r1, genome = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir)
    seq_data.temp_blastn_file = dnaid1315_expected_outcomes/"library_11_1_FKDL202598974-1a-D701-AK1682_HHG5YDSXY_L4_subsample_1.blastn"
    seq_data._find_most_likely_positions(filter_below=0)
    expected_csv = dnaid1315_expected_outcomes/"23-06-22-likely_positions.csv"
    expected_dict = pd.read_csv(expected_csv).set_index('barcode').T.to_dict()
    out_dict = seq_data.positions.set_index('barcode').T.to_dict()
    assert all([v == out_dict[k] for k,v in expected_dict.items()])


def test_merge_colliding_bcs(tn5_structure, map_test_data, tmpdir, dnaid1315_expected_outcomes):
    r1, genome = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir)
    positions = pd.read_csv(dnaid1315_expected_outcomes/"23-06-22-likely_positions.csv")
    seq_data.positions = positions
    seq_data._merge_colliding_barcodes()
    tmp_csv = tmpdir.join('merge_colliding_bcs.csv')
    seq_data.positions.to_csv(tmp_csv, index=False)
    expected_csv = dnaid1315_expected_outcomes / '23-06-22-merge_colliding_bcs.csv'
    assert_files_are_same(tmp_csv, expected_csv)


def test_map(tn5_structure, map_test_data, tmpdir, dnaid1315_expected_outcomes):
    r1, genome = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir, name="23-06-22-library_11_1")
    min_host_bases = 20
    filter_below = 0
    seq_data.map_insertions(min_host_bases, filter_below)
    out_map = tmpdir.join('23-06-22-library_11_1.map.csv')
    expected_map = dnaid1315_expected_outcomes/'23-06-22-library_11_1.map.csv'
    assert_files_are_same(out_map, expected_map)


# def test__find_annotation_overlaps_intersect(sl1344_gff, tmpdir, dnaid1315_expected_outcomes):
#     map_file = dnaid1315_expected_outcomes/'23-06-22-library_11_1.map.csv'
#     annotated_map = AnnotatedMap(map_file, sl1344_gff, "gene", ("ID", "Name", "locus_tag"),
#                                  name='library_11_1', output_dir=tmpdir)
#     annotated_map._find_annotation_overlaps(intersect=True)
#     annotated_map._add_bedtools_results_to_positions(intersect=True)


    # out_intersect = tmpdir.join('library_11_1.bed.intersect.tab')
    # expected_intersect = dnaid1315_expected_outcomes/'25-08-22-library_11_1.bed.intersect.tab' # Cant find this file?
    # assert_files_are_same(out_intersect, expected_intersect)


#
# def test__add_bedtools_results_to_positions(sl1344_gff, tmpdir, dnaid1315_expected_outcomes):
#     map_file = dnaid1315_expected_outcomes/'library_11_1_FKDL202598974-1a-D701-AK1682_HHG5YDSXY_L4_subsample_1.map.csv'
#     annotated_map = AnnotatedMap(map_file, sl1344_gff, "gene", ("ID", "Name", "locus_tag"),
#                                  name='library_11_1', output_dir=tmpdir)
#     annotated_map._find_annotation_overlaps()
#     annotated_map._add_bedtools_results_to_positions(intersect=True)
#     out_map = tmpdir.join('library_11_1.annotated.csv')
#     expected_map = dnaid1315_expected_outcomes/'25-08-22-library_11_1.annotated.csv'
#     assert_files_are_same(out_map, expected_map)
#



# def test__find_annotation_overlaps_closest(sl1344_gff, tmpdir, dnaid1315_expected_outcomes):
#     map_file = dnaid1315_expected_outcomes/'library_11_1_FKDL202598974-1a-D701-AK1682_HHG5YDSXY_L4_subsample_1.map.csv'
#     annotated_map = AnnotatedMap(map_file, sl1344_gff, "gene", ("ID", "Name", "locus_tag"),
#                                  name='library_11_1', output_dir=tmpdir)
#     annotated_map._find_annotation_overlaps(intersect=False)
#     out_intersect = tmpdir.join('library_11_1.bed.intersect.tab')
#     expected_intersect = dnaid1315_expected_outcomes/'25-08-22-library_11_1-closest.bed.intersect.tab'
#     assert_files_are_same(out_intersect, expected_intersect)
#
# def test__add_bedtools_results_to_positions_with_closest(sl1344_gff, tmpdir,dnaid1315_expected_outcomes ):
#     map_file = dnaid1315_expected_outcomes/'library_11_1_FKDL202598974-1a-D701-AK1682_HHG5YDSXY_L4_subsample_1.map.csv'
#     annotated_map = AnnotatedMap(map_file, sl1344_gff, "gene", ("ID", "Name", "locus_tag"),
#                                  name='library_11_1_closest', output_dir=tmpdir)
#     annotated_map._find_annotation_overlaps(intersect=False)
#     annotated_map._add_bedtools_results_to_positions(intersect=False)
#     out_map = tmpdir.join('library_11_1_closest.annotated.csv')
#     expected_map = dnaid1315_expected_outcomes/'25-08-22-library_11_1-closest.annotated.csv'
#     assert_files_are_same(out_map, expected_map)


def test_annotate_intersect(map_test_data, tn5_structure, sl1344_gff, tmpdir,dnaid1315_expected_outcomes):
    map_file = dnaid1315_expected_outcomes / '23-06-22-library_11_1.map.csv'
    feature_type = 'gene'
    identifiers = ('Name', 'locus_tag', 'ID')
    am = AnnotatedMap(map_file=map_file, annotation_file=sl1344_gff,
                                 name='library_11_1',
                                 feature_type=feature_type, identifiers=identifiers, output_dir=tmpdir)
    am.annotate(intersect=True)
    basic = True
    expected = pd.read_json(dnaid1315_expected_outcomes/'annotate_map_expected.json')
    expected['percentile'] = expected['percentile'].round(2)
    to_test = am.annotated_positions[am.annotated_positions.barcode.isin(expected.barcode.values)].copy()

    if basic:
        expected = expected[['barcode', 'chr', 'insertion_site', 'Name', 'locus_tag','distance_to_feature']]
        to_test = to_test[['barcode', 'chr', 'insertion_site', 'Name', 'locus_tag', 'distance_to_feature']]

    expected = expected.set_index('barcode').T.to_dict()
    to_test = to_test.set_index('barcode').T.to_dict()
    #print([k for k,v in expected.items() if v != to_test[k]])
    #print(to_test['CAAGACACTCGATAGAG'])
    #print(expected['CAAGACACTCGATAGAG'])
    assert all([v == to_test[k] for k, v in expected.items()])


def test_annotate_closest(map_test_data, tn5_structure, sl1344_gff, tmpdir,dnaid1315_expected_outcomes):
    map_file = dnaid1315_expected_outcomes / '23-06-22-library_11_1.map.csv'
    feature_type = 'gene'
    identifiers = ('Name','ID')
    am = AnnotatedMap(map_file=map_file, annotation_file=sl1344_gff,
                                 name='library_11_1',
                                 feature_type=feature_type, identifiers=identifiers, output_dir=tmpdir)
    am.annotate(intersect=False)
    expected = pd.read_json(dnaid1315_expected_outcomes/'annotate_map_expected_intersect_false.json')
    expected['distance_to_feature'] = expected['distance_to_feature'].astype(int)
    to_test = am.annotated_positions[am.annotated_positions.barcode.isin(expected.barcode.values)].copy()
    basic=True
    if 'percentile' in to_test.columns:
        assert to_test.percentile.isna().sum() == 100
        to_test = to_test.drop('percentile', axis=1)
    if basic:
        expected = expected[['barcode', 'chr', 'insertion_site', 'Name', 'distance_to_feature']]
        to_test = to_test[['barcode', 'chr', 'insertion_site', 'Name', 'distance_to_feature']]
    expected = expected.set_index('barcode').T.to_dict()
    to_test = to_test.set_index('barcode').T.to_dict()
    print(expected['AAGAGTGCTGGATCCTG'], to_test['AAGAGTGCTGGATCCTG'])
    print ([k for k,v in expected.items() if v != to_test[k]])
    assert all([v == to_test[k] for k, v in expected.items()])



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

