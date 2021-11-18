from mBARq.core import Barcode, BarSeqData


TESTDATA = "./tests/test_files"

# Methods are _parse_structure, extract_barcode_host and editdistance

def get_test_data():
    r1 = f'{TESTDATA}/library_13_1_1.fq'
    r2 = f'{TESTDATA}/library_13_1_2.fq'
    return r1, r2

experiments = ['rbseq', 'wish']

def get_structure(experiment='rbseq'):
    if experiment == 'rbseq':
        return 'GTGTATAAGAGACAG:17:13:before'
    elif experiment == 'wish':
        return 'GGAGGTTCACAATGTGGGAGGTCA:40:0:after'
    else:
        return None


def test__parse_structure():
    # Might completely change how transposon is encoded
    pass


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
