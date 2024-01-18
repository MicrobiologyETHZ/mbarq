from mbarq.core import Barcode, BarSeqData, FastA

# Methods are _parse_structure, extract_barcode_host and editdistance

experiments = ['rbseq', 'wish']

def test__parse_structure():  
    # RBSeq
    barcode = Barcode('B17N13GTGTATAAGAGACAG')
    assert barcode.bc_len == 17
    assert barcode.tn_seq == 'GTGTATAAGAGACAG'
    assert barcode.len_spacer == 13
    assert barcode.bc_before_tn is True

    # WISH
    barcode = Barcode('GGAGGTTCACAATGTGGGAGGTCAB40')
    assert barcode.bc_len == 40
    assert barcode.tn_seq == 'GGAGGTTCACAATGTGGGAGGTCA'
    assert barcode.len_spacer == 0
    assert barcode.bc_before_tn is False

    barcode = Barcode('B20CGTACGCTGCAGGTC')
    assert barcode.bc_len == 20
    assert barcode.tn_seq == 'CGTACGCTGCAGGTC'
    assert barcode.len_spacer == 0
    assert barcode.bc_before_tn is True


def test_extract_barcode_host_rbseq(tn5_structure, map_test_data):
    """
    tn5_structure: 'B17N13GTGTATAAGAGACAG'
    map_test_data: seq_file, genome_file, List[Barcode], mapping_dir

    """

    r1, genome, _, _ = map_test_data
    seq_data = BarSeqData(r1)
    r = list(seq_data.stream_seq_file())[4]
    barcode = Barcode(tn5_structure)
    barcode.extract_barcode_host(r)
    expected_barcode = 'CGTCAGGGCAGCGAACA'
    expected_host = 'AACTCAGTCTAACGCCAAGGGTCTGCTGGGCGCGCTGCGTGATATGCAGGCAAAAGCGAAAGCCGCAGGTCACACGCTGGCGCTCTCCATGAGTATCGGCGGCTG'
    out_barcode, out_host = barcode.bc_seq, barcode.host
    assert out_barcode == expected_barcode
    assert out_host == expected_host


def test_extract_barcode_host_dbauer():

    r1 = FastA('read', "TTGGTGCGCCCTGCAGGGATGTCCACGAGGTCTCTTCTGTACAAAGGGAACTTTGCGTACGCTGCAGGTCGACGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGTTAAAGCCCAAATCAAAATAGTTGGGGCCC")
    barcode = Barcode('B20CGTACGCTGCAGGTC')
    barcode.extract_barcode_host(r1)
    expected_barcode = "TCTGTACAAAGGGAACTTTG"
    expected_host = "GACGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGTTAAAGCCCAAATCAAAATAGTTGGGGCCC"
    out_barcode, out_host = barcode.bc_seq, barcode.host
    assert out_barcode == expected_barcode
    assert out_host == expected_host

    r1 = FastA('read', "GTGTAGCGCCCTGCAGGGATGTCCACGAGGTCTCTCTTATGCTTCACAAATTGAGCGTACGCTGCAGGTCGACGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGTTATCAGTGCCTGCAGAAACCGGGAGATCGGAAGAGCACACG")
    barcode = Barcode('B20N34GACTTATCAGCCAACCTGT')
    barcode.extract_barcode_host(r1)
    expected_barcode = "CTTATGCTTCACAAATTGAG"
    expected_host = "TATCAGTGCCTGCAGAAACCGGGAGATCGGAAGAGCACACG"

    out_barcode, out_host = barcode.bc_seq, barcode.host
    assert out_barcode == expected_barcode
    assert out_host == expected_host

# todo add tests for WISH barcode structure

# todo add tests for mariner barcode structure
