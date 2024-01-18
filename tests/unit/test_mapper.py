from mbarq.mapper import Mapper, AnnotatedMap
from mbarq.core import Barcode
from collections import Counter
from pathlib import Path
import pandas as pd
import pickle
from test_utils import assert_files_are_same
import numpy as np


def test_extract_barcodes_tn5(tn5_structure, map_test_data, tmpdir):
    """Check ...."""
    r1, genome, map_barcodes, _ = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir)
    seq_data.extract_barcodes()
    barcodes = [(bc.bc_seq, bc.host) for bc in seq_data.barcodes]
    expected_barcodes = [(bc.bc_seq, bc.host) for bc in map_barcodes]
    assert Counter(barcodes) == Counter(expected_barcodes)


def test__dereplicate_barcodes(tn5_structure, map_test_data, tmpdir):
    r1, genome, map_barcodes, _ = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir)
    seq_data.barcodes = map_barcodes
    barcodes = [barcode for barcode in seq_data._dereplicate_barcodes()]
    assert len(barcodes) == 5
    assert barcodes[0].count == 2


def test__write_barcodes_to_fasta(tn5_structure, map_test_data, tmpdir):
    r1, genome, map_barcodes, mapping_dir = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome,
                      name='mapping_library', output_dir=tmpdir)
    seq_data.barcodes = map_barcodes
    seq_data._write_barcodes_to_fasta()
    out_fasta = tmpdir.join("mapping_library.fasta")
    expected_fasta_file = mapping_dir / "mapping_library.fasta"
    assert_files_are_same(out_fasta, expected_fasta_file)


def test__blast_barcode_host_genome(tn5_structure, map_test_data, tmpdir):
    r1, genome, map_barcodes, mapping_dir = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir,
                      name='mapping_library')
    seq_data.barcodes = map_barcodes
    seq_data._write_barcodes_to_fasta()
    seq_data._blast_barcode_host()
    out_blastn_file = tmpdir.join("mapping_library.blastn")
    expected_blastn_file = mapping_dir / "mapping_library.blastn"
    assert_files_are_same(out_blastn_file, expected_blastn_file)

    # todo clean up db index


def test__blast_barcode_host_genome_gzip(tn5_structure, map_test_data, tmpdir):
    r1, genome, map_barcodes, mapping_dir = map_test_data
    genome = mapping_dir/"test_genome_gz.fna.gz"
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir,
                      name='mapping_library')
    seq_data.barcodes = map_barcodes
    seq_data._write_barcodes_to_fasta()
    seq_data._blast_barcode_host()
    out_blastn_file = tmpdir.join("mapping_library.blastn")
    expected_blastn_file = mapping_dir / "mapping_library.blastn"
    assert_files_are_same(out_blastn_file, expected_blastn_file)
    assert Path(genome).is_file()

# todo add test for multimap removal


def test__find_most_likely_positions(tn5_structure, map_test_data, tmpdir):
    r1, genome, mapping_barcodes, mapping_dir = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir)
    seq_data.temp_blastn_file = mapping_dir / "mapping_library.blastn"
    seq_data._find_most_likely_positions(filter_below=0)
    expected_dict = {
        'ATGAAAAGTATGAATCC': {'sseqid': 'FQ312003.1', 'sstrand': 'minus', 'insertion_site': 16097, 'total_count': 3,
                              'prop_read_per_position': 1.0},
        'CGTCAGGGCAGCGAACA': {'sseqid': 'FQ312003.1', 'sstrand': 'plus', 'insertion_site': 18337, 'total_count': 2,
                              'prop_read_per_position': 1.0},
        'TTGGGATCCCACCATTT': {'sseqid': 'FQ312003.1', 'sstrand': 'plus', 'insertion_site': 17682, 'total_count': 1,
                              'prop_read_per_position': 1.0}}
    out_dict = seq_data.positions.set_index('barcode').T.to_dict()
    assert all([v == out_dict[k] for k, v in expected_dict.items()])


# def test_merge_colliding_bcs(tn5_structure, map_test_data, tmpdir):
#     r1, genome, _, _ = map_test_data
#     seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir)
#     positions = """AAAACCTCCCTACCCAT,HE654726.1,minus,7626,1,1.0
#         AAAACCTCCCTGCCCAT,HE654726.1,minus,7626,87,1.0
#         AAAACCTCCCTGTCCAT,HE654726.1,minus,7626,1,1.0
#         AGGACCTCCCTGGGGAT,HE654726.1,minus,7626,1,1.0
#     """.split()
#     positions = (pd.DataFrame([p.split(",") for p in positions],
#                               columns=['barcode', 'sseqid', 'sstrand', 'insertion_site',
#                                        'total_count', 'prop_read_per_position'])
#                  .astype({'insertion_site': int,
#                           'total_count': int,
#                           'prop_read_per_position': float}))
#     seq_data.positions = positions
#     seq_data._merge_colliding_barcodes()
#     out_dict = seq_data.positions.set_index('barcode').to_dict()
#     expected_dict = {'total_count': {'AAAACCTCCCTGCCCAT': 89, 'AGGACCTCCCTGGGGAT': 1},
#                      'insertion_site': {'AAAACCTCCCTGCCCAT': 7626, 'AGGACCTCCCTGGGGAT': 7626},
#                      'sseqid': {'AAAACCTCCCTGCCCAT': 'HE654726.1', 'AGGACCTCCCTGGGGAT': 'HE654726.1'},
#                      'sstrand': {'AAAACCTCCCTGCCCAT': 'minus', 'AGGACCTCCCTGGGGAT': 'minus'}}
#     assert all([v == out_dict[k] for k, v in expected_dict.items()])
#     assert (len(out_dict) == len(expected_dict))


def test_merge_colliding_bcs(tn5_structure, map_test_data, tmpdir):
    r1, genome, _, _ = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir)

    # Test1
    print('Test 1')
    seq_data.positions = pd.DataFrame({'barcode': {0: 'AAAACCTCCCTACCCAT',
                                                   1: 'AAAACCTCCCTGCCCAT',
                                                   2: 'AAAACCTCCCTGTCCAT',
                                                   3: 'AGGACCTCCCTGGGGAT'},
                                       'sseqid': {0: 'HE654726.1',
                                                  1: 'HE654726.1',
                                                  2: 'HE654726.1',
                                                  3: 'HE654726.1'},
                                       'sstrand': {0: 'minus', 1: 'minus', 2: 'minus', 3: 'minus'},
                                       'insertion_site': {0: 7626, 1: 7626, 2: 7626, 3: 7626},
                                       'total_count': {0: 1, 1: 87, 2: 1, 3: 1},
                                       'prop_read_per_position': {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0}})
    seq_data._merge_colliding_barcodes()
    expected_positions = pd.DataFrame({'barcode': {0: 'AAAACCTCCCTGCCCAT', 3: 'AGGACCTCCCTGGGGAT'},
                                       'total_count': {0: 89, 3: 1},
                                       'insertion_site': {0: 7626, 3: 7626},
                                       'sseqid': {0: 'HE654726.1', 3: 'HE654726.1'},
                                       'sstrand': {0: 'minus', 3: 'minus'}})
    print(seq_data.positions)
    print(expected_positions)
    assert expected_positions.equals(seq_data.positions)

    # Test 2
    print('Test 2')
    seq_data.positions = pd.DataFrame({'barcode': {16182: 'ACCACAACCACGCTCAG',
                                                   16184: 'ACCCCAACCACGCTCAG',
                                                   16183: 'CGCAGCTATCCAACCCA'},
                                       'total_count': {16182: 4, 16184: 92, 16183: 68},
                                       'insertion_site': {16182: 1814830, 16184: 1814830, 16183: 1814830},
                                       'sseqid': {16182: 1, 16184: 1, 16183: 1},
                                       'sstrand': {16182: 'minus', 16184: 'minus', 16183: 'minus'}})
    expected_positions = pd.DataFrame({'barcode': {16182: 'ACCCCAACCACGCTCAG', 16183: 'CGCAGCTATCCAACCCA'},
                                       'total_count': {16182: 96, 16183: 68},
                                       'insertion_site': {16182: 1814830, 16183: 1814830},
                                       'sseqid': {16182: 1, 16183: 1},
                                       'sstrand': {16182: 'minus', 16183: 'minus'}
                                       })
    seq_data._merge_colliding_barcodes()

    assert expected_positions.equals(seq_data.positions)

    # Test 3
    print('Test 3')
    seq_data.positions = pd.DataFrame({'barcode': {16797: 'ACAGGCACACCCAACAC',
                                                   16795: 'ACAGGTACACCCAACAC',
                                                   16796: 'CACCACACTCCAAAACC'},
                                       'total_count': {16797: 4, 16795: 138, 16796: 18},
                                       'insertion_site': {16797: 2206977, 16795: 2206977, 16796: 2206977},
                                       'sseqid': {16797: 1, 16795: 1, 16796: 1},
                                       'sstrand': {16797: 'plus', 16795: 'plus', 16796: 'plus'}})
    expected_positions = pd.DataFrame({'barcode': {16797: 'ACAGGTACACCCAACAC', 16796: 'CACCACACTCCAAAACC'},
                                       'total_count': {16797: 142, 16796: 18},
                                       'insertion_site': {16797: 2206977, 16796: 2206977},
                                       'sseqid': {16797: 1, 16796: 1},
                                       'sstrand': {16797: 'plus',  16796: 'plus'}})
    seq_data._merge_colliding_barcodes()
    print(seq_data.positions)
    print(expected_positions)
    assert expected_positions.equals(seq_data.positions)

    # Test 4
    print('Test 4')
    seq_data.positions = pd.DataFrame({'barcode': {27210: 'CAACTAAGCCAAAAAGG',
                                                   27212: 'CAACTAAGCCAAACAGG',
                                                   27213: 'CTACATTGAGCCAAGCC',
                                                   27211: 'CTACATTGCGCCAAGCC',
                                                   27214: 'GAGACTCAGCTGTTGCG',
                                                   27215: 'GCGATCAGCGAAAAAAG'},
                                       'total_count': {27210: 1,
                                                       27212: 130,
                                                       27213: 1,
                                                       27211: 59,
                                                       27214: 124,
                                                       27215: 18},
                                       'insertion_site': {27210: 107077,
                                                          27212: 107077,
                                                          27213: 107077,
                                                          27211: 107077,
                                                          27214: 107077,
                                                          27215: 107077},
                                       'sseqid': {27210: 2, 27212: 2, 27213: 2, 27211: 2, 27214: 2, 27215: 2},
                                       'sstrand': {27210: 'minus',
                                                   27212: 'minus',
                                                   27213: 'minus',
                                                   27211: 'minus',
                                                   27214: 'minus',
                                                   27215: 'minus'}})
    expected_positions = pd.DataFrame({'barcode': {27210: 'CAACTAAGCCAAACAGG',
                                                   27213: 'CTACATTGCGCCAAGCC',
                                                   27214: 'GAGACTCAGCTGTTGCG',
                                                   27215: 'GCGATCAGCGAAAAAAG'},
                                       'total_count': {27210: 131,
                                                       27213: 60,
                                                       27214: 124,
                                                       27215: 18},
                                       'insertion_site': {27210: 107077,
                                                          27213: 107077,
                                                          27214: 107077,
                                                          27215: 107077},
                                       'sseqid': {27210: 2, 27213: 2, 27214: 2, 27215: 2},
                                       'sstrand': {27210: 'minus',
                                                   27213: 'minus',
                                                   27214: 'minus',
                                                   27215: 'minus',
                                                   }})
    seq_data._merge_colliding_barcodes()
    assert expected_positions.equals(seq_data.positions)

    # Test 5
    print('Test 5')
    seq_data.positions = pd.DataFrame({'barcode': {13908: 'CAACAAGAACACACAAG',
                                                   13909: 'AAAGGACATACAAAGAC',
                                                   13910: 'TCAGCCAACGGAACAAG',
                                                   13911: 'CAACAAGAACTCACAAG',
                                                   13912: 'AAAGGACATACAAAGAA'},
                                       'total_count': {13908: 350, 13909: 1, 13910: 1, 13911: 3, 13912: 78},
                                       'insertion_site': {13908: 881783,
                                                          13909: 881783,
                                                          13910: 881783,
                                                          13911: 881783,
                                                          13912: 881783},
                                       'sseqid': {13908: 1, 13909: 1, 13910: 1, 13911: 1, 13912: 1},
                                       'sstrand': {13908: 'minus',
                                                   13909: 'minus',
                                                   13910: 'minus',
                                                   13911: 'minus',
                                                   13912: 'minus'}})
    expected_positions = pd.DataFrame({'barcode': {13912: 'AAAGGACATACAAAGAA',
                                                   13908: 'CAACAAGAACACACAAG',
                                                   13910: 'TCAGCCAACGGAACAAG'},
                                       'total_count': {13912: 79, 13908: 353, 13910: 1},
                                       'insertion_site': {13912: 881783, 13908: 881783, 13910: 881783},
                                       'sseqid': {13912: 1, 13908: 1, 13910: 1},
                                       'sstrand': {13912: 'minus', 13908: 'minus', 13910: 'minus'}})
    seq_data._merge_colliding_barcodes()
    assert expected_positions.equals(seq_data.positions)

    # Test 6
    print('Test 6')
    seq_data.positions = pd.DataFrame({'barcode': {27095: 'CAATATTTAAGCAGCAC',
                                                   27096: 'ACACCGCAAACAAAGGA',
                                                   27097: 'CAATATTTAAGCCGCAC'},
                                       'total_count': {27095: 1, 27096: 92, 27097: 69},
                                       'insertion_site': {27095: 103567, 27096: 103567, 27097: 103567},
                                       'sseqid': {27095: 2, 27096: 2, 27097: 2},
                                       'sstrand': {27095: 'minus', 27096: 'minus', 27097: 'minus'}})
    expected_positions = pd.DataFrame({'barcode': {27096: 'ACACCGCAAACAAAGGA', 27095: 'CAATATTTAAGCCGCAC'},
                                       'total_count': {27096: 92, 27095: 70},
                                       'insertion_site': {27096: 103567, 27095: 103567},
                                       'sseqid': {27096: 2, 27095: 2},
                                       'sstrand': {27096: 'minus', 27095: 'minus'}})
    seq_data._merge_colliding_barcodes()
    assert expected_positions.equals(seq_data.positions)

    # Test 7
    print('Test 7')
    seq_data.positions = pd.DataFrame({'barcode': {26600: 'CTCAACCACCAACTGTG',
                                                   26601: 'TCATCGAGCAGCCAAAC',
                                                   26602: 'CTCAACCACCAACGGTG',
                                                   26603: 'AAAGAGCCAAAACCATA'},
                                       'total_count': {26600: 2, 26601: 151, 26602: 81, 26603: 35},
                                       'insertion_site': {26600: 23803, 26601: 23803, 26602: 23803, 26603: 23803},
                                       'sseqid': {26600: 2, 26601: 2, 26602: 2, 26603: 2},
                                       'sstrand': {26600: 'minus', 26601: 'minus', 26602: 'minus', 26603: 'minus'}})
    expected_positions = pd.DataFrame({'barcode': {26603: 'AAAGAGCCAAAACCATA',
                                                   26602: 'CTCAACCACCAACGGTG',
                                                   26601: 'TCATCGAGCAGCCAAAC'},
                                       'total_count': {26603: 35, 26602: 83, 26601: 151},
                                       'insertion_site': {26603: 23803, 26602: 23803, 26601: 23803},
                                       'sseqid': {26603: 2, 26602: 2, 26601: 2},
                                       'sstrand': {26603: 'minus', 26602: 'minus', 26601: 'minus'}})
    seq_data._merge_colliding_barcodes()
    assert expected_positions.equals(seq_data.positions)


def test_map(tn5_structure, map_test_data, tmpdir):
    r1, genome, _, _ = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir,
                      name="test_map")
    min_host_bases = 20
    filter_below = 0
    seq_data.map_insertions(min_host_bases, filter_below)
    out_dict = seq_data.positions.set_index('barcode').to_dict()

    expected_dict = {
        'abundance_in_mapping_library': {'ATGAAAAGTATGAATCC': 3, 'TTGGGATCCCACCATTT': 1, 'CGTCAGGGCAGCGAACA': 2},
        'insertion_site': {'ATGAAAAGTATGAATCC': 16097, 'TTGGGATCCCACCATTT': 17682, 'CGTCAGGGCAGCGAACA': 18337},
        'chr': {'ATGAAAAGTATGAATCC': 'FQ312003.1', 'TTGGGATCCCACCATTT': 'FQ312003.1',
                'CGTCAGGGCAGCGAACA': 'FQ312003.1'},
        'strand': {'ATGAAAAGTATGAATCC': 'minus', 'TTGGGATCCCACCATTT': 'plus', 'CGTCAGGGCAGCGAACA': 'plus'}}

    assert all([v == out_dict[k] for k, v in expected_dict.items()])
    assert (len(out_dict) == len(expected_dict))


# ANNOTATION


def test__find_annotation_overlaps(test_gff, map_test_data, tmpdir):
    _, _, _, mapping_dir = map_test_data
    map_file = mapping_dir / 'test_map.map.csv'
    annotated_map = AnnotatedMap(map_file, test_gff, "gene", ("ID", "Name", "locus_tag"),
                                 name='test_library', output_dir=tmpdir)
    annotated_map._find_annotation_overlaps(intersect=False)
    out_dict = annotated_map.annotated_positions.set_index(
        'barcode').fillna('nan').to_dict()
    expected_dict = {'chr': {'ATGAAAAGTATGAATCC': 'FQ312003.1', 'TTGGGATCCCACCATTT': 'FQ312003.1',
                             'CGTCAGGGCAGCGAACA': 'FQ312003.1'},
                     'insertion_site': {'ATGAAAAGTATGAATCC': 16097, 'TTGGGATCCCACCATTT': 17682,
                                        'CGTCAGGGCAGCGAACA': 18337},
                     'abundance_in_mapping_library': {'ATGAAAAGTATGAATCC': 3, 'TTGGGATCCCACCATTT': 1,
                                                      'CGTCAGGGCAGCGAACA': 2},
                     'gene_start': {'ATGAAAAGTATGAATCC': 16088, 'TTGGGATCCCACCATTT': 17043, 'CGTCAGGGCAGCGAACA': 18083},
                     'gene_end': {'ATGAAAAGTATGAATCC': 16432, 'TTGGGATCCCACCATTT': 17486, 'CGTCAGGGCAGCGAACA': 19966},
                     'gene_strand': {'ATGAAAAGTATGAATCC': '+', 'TTGGGATCCCACCATTT': '-', 'CGTCAGGGCAGCGAACA': '+'},
                     'ID': {'ATGAAAAGTATGAATCC': 'gene-SL1344_0015', 'TTGGGATCCCACCATTT': 'gene-SL1344_0017',
                            'CGTCAGGGCAGCGAACA': 'gene-SL1344_0018'},
                     'Name': {'ATGAAAAGTATGAATCC': 'SL1344_0015', 'TTGGGATCCCACCATTT': 'SL1344_0017',
                              'CGTCAGGGCAGCGAACA': 'chiA'},
                     'locus_tag': {'ATGAAAAGTATGAATCC': 'SL1344_0015', 'TTGGGATCCCACCATTT': 'SL1344_0017',
                                   'CGTCAGGGCAGCGAACA': 'SL1344_0018'},
                     'distance_to_feature': {'ATGAAAAGTATGAATCC': 0, 'TTGGGATCCCACCATTT': -196, 'CGTCAGGGCAGCGAACA': 0},
                     'percentile': {'ATGAAAAGTATGAATCC': 0.03, 'TTGGGATCCCACCATTT': 'nan', 'CGTCAGGGCAGCGAACA': 0.13}}
    assert all([v == out_dict[k] for k, v in expected_dict.items()])
    assert (len(out_dict) == len(expected_dict))

    annotated_map._find_annotation_overlaps(intersect=True)
    expected_dict = {'chr': {'ATGAAAAGTATGAATCC': 'FQ312003.1', 'CGTCAGGGCAGCGAACA': 'FQ312003.1'},
                     'insertion_site': {'ATGAAAAGTATGAATCC': 16097, 'CGTCAGGGCAGCGAACA': 18337},
                     'abundance_in_mapping_library': {'ATGAAAAGTATGAATCC': 3, 'CGTCAGGGCAGCGAACA': 2},
                     'gene_start': {'ATGAAAAGTATGAATCC': 16088, 'CGTCAGGGCAGCGAACA': 18083},
                     'gene_end': {'ATGAAAAGTATGAATCC': 16432, 'CGTCAGGGCAGCGAACA': 19966},
                     'gene_strand': {'ATGAAAAGTATGAATCC': '+', 'CGTCAGGGCAGCGAACA': '+'},
                     'ID': {'ATGAAAAGTATGAATCC': 'gene-SL1344_0015', 'CGTCAGGGCAGCGAACA': 'gene-SL1344_0018'},
                     'Name': {'ATGAAAAGTATGAATCC': 'SL1344_0015', 'CGTCAGGGCAGCGAACA': 'chiA'},
                     'locus_tag': {'ATGAAAAGTATGAATCC': 'SL1344_0015', 'CGTCAGGGCAGCGAACA': 'SL1344_0018'},
                     'distance_to_feature': {'ATGAAAAGTATGAATCC': 0, 'CGTCAGGGCAGCGAACA': 0},
                     'percentile': {'ATGAAAAGTATGAATCC': 0.03, 'CGTCAGGGCAGCGAACA': 0.13}}
    out_dict = annotated_map.annotated_positions.set_index('barcode').to_dict()
    assert all([v == out_dict[k] for k, v in expected_dict.items()])
    assert (len(out_dict) == len(expected_dict))


def test_annotate_intersect(map_test_data, test_gff, tmpdir):
    _, _, _, mapping_dir = map_test_data
    map_file = mapping_dir / '23-06-22-library_11_1.map.csv'
    feature_type = 'gene'
    identifiers = ('Name', 'locus_tag', 'ID')
    am = AnnotatedMap(map_file=map_file, annotation_file=test_gff,
                      name='test_library',
                      feature_type=feature_type, identifiers=identifiers, output_dir=tmpdir)
    am.annotate(intersect=True)
    expected_dict = {'chr': {'ATGAAAAGTATGAATCC': 'FQ312003.1', 'CGTCAGGGCAGCGAACA': 'FQ312003.1',
                             'TGTATTGCCCAATAATG': 'FQ312003.1'},
                     'insertion_site': {'ATGAAAAGTATGAATCC': 16097, 'CGTCAGGGCAGCGAACA': 18337,
                                        'TGTATTGCCCAATAATG': 19222},
                     'abundance_in_mapping_library': {'ATGAAAAGTATGAATCC': 4, 'CGTCAGGGCAGCGAACA': 4,
                                                      'TGTATTGCCCAATAATG': 4},
                     'gene_start': {'ATGAAAAGTATGAATCC': 16088, 'CGTCAGGGCAGCGAACA': 18083, 'TGTATTGCCCAATAATG': 18083},
                     'gene_end': {'ATGAAAAGTATGAATCC': 16432, 'CGTCAGGGCAGCGAACA': 19966, 'TGTATTGCCCAATAATG': 19966},
                     'gene_strand': {'ATGAAAAGTATGAATCC': '+', 'CGTCAGGGCAGCGAACA': '+', 'TGTATTGCCCAATAATG': '+'},
                     'Name': {'ATGAAAAGTATGAATCC': 'SL1344_0015', 'CGTCAGGGCAGCGAACA': 'chiA',
                              'TGTATTGCCCAATAATG': 'chiA'},
                     'locus_tag': {'ATGAAAAGTATGAATCC': 'SL1344_0015', 'CGTCAGGGCAGCGAACA': 'SL1344_0018',
                                   'TGTATTGCCCAATAATG': 'SL1344_0018'},
                     'ID': {'ATGAAAAGTATGAATCC': 'gene-SL1344_0015', 'CGTCAGGGCAGCGAACA': 'gene-SL1344_0018',
                            'TGTATTGCCCAATAATG': 'gene-SL1344_0018'},
                     'distance_to_feature': {'ATGAAAAGTATGAATCC': 0, 'CGTCAGGGCAGCGAACA': 0, 'TGTATTGCCCAATAATG': 0},
                     'percentile': {'ATGAAAAGTATGAATCC': 0.03, 'CGTCAGGGCAGCGAACA': 0.13, 'TGTATTGCCCAATAATG': 0.6}}

    out_dict = am.annotated_positions.set_index('barcode').to_dict()
    assert all([v == out_dict[k] for k, v in expected_dict.items()])
    assert (len(out_dict) == len(expected_dict))


def test_annotate_closest(map_test_data, tn5_structure, test_gff, tmpdir):
    _, _, _, mapping_dir = map_test_data
    map_file = mapping_dir / '23-06-22-library_11_1.map.csv'
    feature_type = 'gene'
    identifiers = ('Name', 'ID')
    am = AnnotatedMap(map_file=map_file, annotation_file=test_gff,
                      name='test_library',
                      feature_type=feature_type, identifiers=identifiers, output_dir=tmpdir)
    am.annotate(intersect=False)

    out_dict = am.annotated_positions
    assert out_dict.Name.isna().sum() == 181
    out_dict = out_dict.dropna(subset=['Name']).fillna('nan')
    out_dict = out_dict[abs(out_dict.distance_to_feature)
                        < 3000].set_index('barcode').to_dict()
    expected_dict = {
        'chr': {'CGTGACAAGCCACTCGG': 'FQ312003.1', 'CCTGTCGCGATAGCTTG': 'FQ312003.1', 'ATGAAAAGTATGAATCC': 'FQ312003.1',
                'TTGGGATCCCACCATTT': 'FQ312003.1', 'CGTCAGGGCAGCGAACA': 'FQ312003.1', 'TGTATTGCCCAATAATG': 'FQ312003.1',
                'TGTGCGGTGGCTATCAC': 'FQ312003.1'},
        'insertion_site': {'CGTGACAAGCCACTCGG': 15210, 'CCTGTCGCGATAGCTTG': 15213, 'ATGAAAAGTATGAATCC': 16097,
                           'TTGGGATCCCACCATTT': 17682, 'CGTCAGGGCAGCGAACA': 18337, 'TGTATTGCCCAATAATG': 19222,
                           'TGTGCGGTGGCTATCAC': 21143},
        'abundance_in_mapping_library': {'CGTGACAAGCCACTCGG': 8, 'CCTGTCGCGATAGCTTG': 3, 'ATGAAAAGTATGAATCC': 4,
                                         'TTGGGATCCCACCATTT': 1, 'CGTCAGGGCAGCGAACA': 4, 'TGTATTGCCCAATAATG': 4,
                                         'TGTGCGGTGGCTATCAC': 1},
        'gene_start': {'CGTGACAAGCCACTCGG': 16088, 'CCTGTCGCGATAGCTTG': 16088, 'ATGAAAAGTATGAATCC': 16088,
                       'TTGGGATCCCACCATTT': 17043, 'CGTCAGGGCAGCGAACA': 18083, 'TGTATTGCCCAATAATG': 18083,
                       'TGTGCGGTGGCTATCAC': 18083},
        'gene_end': {'CGTGACAAGCCACTCGG': 16432, 'CCTGTCGCGATAGCTTG': 16432, 'ATGAAAAGTATGAATCC': 16432,
                     'TTGGGATCCCACCATTT': 17486, 'CGTCAGGGCAGCGAACA': 19966, 'TGTATTGCCCAATAATG': 19966,
                     'TGTGCGGTGGCTATCAC': 19966},
        'gene_strand': {'CGTGACAAGCCACTCGG': '+', 'CCTGTCGCGATAGCTTG': '+', 'ATGAAAAGTATGAATCC': '+',
                        'TTGGGATCCCACCATTT': '-', 'CGTCAGGGCAGCGAACA': '+', 'TGTATTGCCCAATAATG': '+',
                        'TGTGCGGTGGCTATCAC': '+'},
        'Name': {'CGTGACAAGCCACTCGG': 'SL1344_0015', 'CCTGTCGCGATAGCTTG': 'SL1344_0015',
                 'ATGAAAAGTATGAATCC': 'SL1344_0015', 'TTGGGATCCCACCATTT': 'SL1344_0017', 'CGTCAGGGCAGCGAACA': 'chiA',
                 'TGTATTGCCCAATAATG': 'chiA', 'TGTGCGGTGGCTATCAC': 'chiA'},
        'ID': {'CGTGACAAGCCACTCGG': 'gene-SL1344_0015', 'CCTGTCGCGATAGCTTG': 'gene-SL1344_0015',
               'ATGAAAAGTATGAATCC': 'gene-SL1344_0015', 'TTGGGATCCCACCATTT': 'gene-SL1344_0017',
               'CGTCAGGGCAGCGAACA': 'gene-SL1344_0018', 'TGTATTGCCCAATAATG': 'gene-SL1344_0018',
               'TGTGCGGTGGCTATCAC': 'gene-SL1344_0018'},
        'distance_to_feature': {'CGTGACAAGCCACTCGG': -878, 'CCTGTCGCGATAGCTTG': -875, 'ATGAAAAGTATGAATCC': 0,
                                'TTGGGATCCCACCATTT': -196, 'CGTCAGGGCAGCGAACA': 0, 'TGTATTGCCCAATAATG': 0,
                                'TGTGCGGTGGCTATCAC': 1177},
        'percentile': {'CGTGACAAGCCACTCGG': 'nan', 'CCTGTCGCGATAGCTTG': 'nan', 'ATGAAAAGTATGAATCC': 0.03,
                       'TTGGGATCCCACCATTT': 'nan', 'CGTCAGGGCAGCGAACA': 0.13, 'TGTATTGCCCAATAATG': 0.6,
                       'TGTGCGGTGGCTATCAC': 'nan'}}
    assert (len(out_dict) == len(expected_dict))
    assert all([v == out_dict[k] for k, v in expected_dict.items()])


def test_annotate_closest_numeric_chr(map_test_data, tmpdir):
    _, _, _, mapping_dir = map_test_data
    map_file = mapping_dir / "numeric_chr_map.csv"
    gff_file = mapping_dir / 'numeric_chr.gff'
    feature_type = 'CDS'
    identifiers = ('ID', 'locus_tag')
    am = AnnotatedMap(map_file=map_file, annotation_file=gff_file,
                      name='test_numeric_chr',
                      feature_type=feature_type, identifiers=identifiers, output_dir=tmpdir)
    am._validate_annotations()
    am._find_annotation_overlaps(intersect=False)
    assert am.annotated_positions['ID'].isna().sum() == 0
    assert am.annotated_positions['locus_tag'].isna().sum() == 0
