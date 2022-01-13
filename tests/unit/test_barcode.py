from mbarq.core import Barcode, BarSeqData


TESTDATA = "./tests/test_files"

# Methods are _parse_structure, extract_barcode_host and editdistance

def get_test_data():
    r1 = f'{TESTDATA}/library_13_1_1.fq'
    r2 = f'{TESTDATA}/library_13_1_2.fq'
    return r1, r2

experiments = ['rbseq', 'wish']

def get_structure(experiment='rbseq'):
    if experiment == 'rbseq':
        #return 'GTGTATAAGAGACAG:17:13:before'
        return 'B17N13GTGTATAAGAGACAG'
    elif experiment == 'wish':
        return 'GGAGGTTCACAATGTGGGAGGTCAB40'
    else:
        return None


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


def test_extract_barcode_host_rbseq():
    structure = get_structure('rbseq')
    r1, r2 = get_test_data()
    seq_data = BarSeqData(r1)
    r = list(seq_data.stream_seq_file())[1000]
    barcode = Barcode(structure)
    barcode.extract_barcode_host(r)
    expected_barcode = 'GACGGCTATACTCAAAG'
    expected_host = 'GTCCTGAACTCGCGACGCAAATAAACGGTATCTTGAGGTATCTAATAACGAGAAATGCTTTAAAATATTAATTTCCGGGTAAGCATCATCATCATAGAAAAATAC'
    out_barcode, out_host = barcode.bc_seq, barcode.host
    assert out_barcode == expected_barcode
    assert out_host == expected_host


# todo add tests for WISH barcode structure
# Leave till final decision made on how to encode a transposon
