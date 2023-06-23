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
    seq_data = Mapper(r1, tn5_structure, genome=genome, name='mapping_library', output_dir=tmpdir)
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


def test_merge_colliding_bcs(tn5_structure, map_test_data, tmpdir):
    r1, genome, _, _ = map_test_data
    seq_data = Mapper(r1, tn5_structure, genome=genome, output_dir=tmpdir)
    positions = """AAAACCTCCCTACCCAT,HE654726.1,minus,7626,1,1.0
        AAAACCTCCCTGCCCAT,HE654726.1,minus,7626,87,1.0
        AAAACCTCCCTGTCCAT,HE654726.1,minus,7626,1,1.0
        AGGACCTCCCTGGGGAT,HE654726.1,minus,7626,1,1.0
    """.split()
    positions = (pd.DataFrame([p.split(",") for p in positions],
                              columns=['barcode', 'sseqid', 'sstrand', 'insertion_site',
                                       'total_count', 'prop_read_per_position'])
                 .astype({'insertion_site': int,
                          'total_count': int,
                          'prop_read_per_position': float}))
    seq_data.positions = positions
    seq_data._merge_colliding_barcodes()
    out_dict = seq_data.positions.set_index('barcode').to_dict()
    expected_dict = {'total_count': {'AAAACCTCCCTGCCCAT': 89, 'AGGACCTCCCTGGGGAT': 1},
                     'insertion_site': {'AAAACCTCCCTGCCCAT': 7626, 'AGGACCTCCCTGGGGAT': 7626},
                     'sseqid': {'AAAACCTCCCTGCCCAT': 'HE654726.1', 'AGGACCTCCCTGGGGAT': 'HE654726.1'},
                     'sstrand': {'AAAACCTCCCTGCCCAT': 'minus', 'AGGACCTCCCTGGGGAT': 'minus'}}
    assert all([v == out_dict[k] for k, v in expected_dict.items()])
    assert (len(out_dict) == len(expected_dict))


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
    out_dict = annotated_map.annotated_positions.set_index('barcode').fillna('nan').to_dict()
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
    out_dict = out_dict[abs(out_dict.distance_to_feature) < 3000].set_index('barcode').to_dict()
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


def test_annotate_closest_numeric_chr(map_test_data):
    _, _, _, mapping_dir = map_test_data
    map_file = mapping_dir / "numeric_chr_map.csv"
    gff_file = mapping_dir / 'numeric_chr.gff'
    feature_type = 'CDS'
    identifiers = ('ID', 'locus_tag')
    am = AnnotatedMap(map_file=map_file, annotation_file=gff_file,
                      name='test_numeric_chr',
                      feature_type=feature_type, identifiers=identifiers, output_dir=".")
    am._validate_annotations()
    am._find_annotation_overlaps(intersect=False)
    assert am.annotated_positions['ID'].isna().sum() == 0
    assert am.annotated_positions['locus_tag'].isna().sum() == 0
