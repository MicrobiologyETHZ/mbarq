from mbarq.core import Barcode, BarSeqData

# Methods are _parse_structure, extract_barcode_host and editdistance

experiments = ['rbseq', 'wish']

def test__parse_structure():  # todo !!!!
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


def test_extract_barcode_host_rbseq(tn5_structure, map_test_data):
    r1, genome = map_test_data
    seq_data = BarSeqData(r1)
    r = list(seq_data.stream_seq_file())[1000]
    barcode = Barcode(tn5_structure)
    barcode.extract_barcode_host(r)
    expected_barcode = 'GCGGTGACAGGATCGGA'
    expected_host = 'AGCAACAACTCTGCTATGAGATTTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATTACTCGATCGCGTATGCCGTCTTCTGCTTGAAAAGGGGGGGGGGG'
    out_barcode, out_host = barcode.bc_seq, barcode.host
    assert out_barcode == expected_barcode
    assert out_host == expected_host


# todo add tests for WISH barcode structure

# todo add tests for mariner barcode structure
