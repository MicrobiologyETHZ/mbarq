from collections import Counter
from Bio import SeqIO
from tnseq import sequence
from tnseq import tnseq
from tnseq import extract_barcodes
from tnseq import quantify_barcodes
import json
import subprocess
import shlex




def test_wishtag_extraction():
    in_r1_file = "/nfs/home/ansintsova/wishtags/MixA_S1_L001_R1_001.fastq"
    r = sequence.stream_fa(in_r1_file)
    tnSeq = 'GGAGGTTCACAATGTGGGAGGTCA'
    #tnSeq = 'TGACCTCCCACATTGTGAACCTCC'
    barcodes = []
    hosts = []
    c = 0
    for seq in r:
        if tnSeq in seq.sequence:
            barcode, host = extract_barcodes.extract_barcode_host(seq, tnSeq, bc2tp2=0, bcLen=40, before=False)
            barcodes.append(barcode)
            hosts.append(host)
        else:
            c += 1
    print(c)
    print(len(barcodes))
    print(barcodes[0:10])
    print(len(barcodes[0]))

#test_wishtag_extraction()



# def test_quantify_load_fq_barcodes():
#     in_fq_file = "/nfs/home/ansintsova/wishtags/MixA_S1_L001_R1_001.fastq"
#     tnSeq = 'GGAGGTTCACAATGTGGGAGGTCA'
#
#     cnter = quantify_barcodes.quantify_load_fq_barcodes(in_fq_file,   tp2=tnSeq, bc2tp2=0, bcLen=40, before=False)
#     print(cnter.most_common()[0:10])
#
# #test_quantify_load_fq_barcodes()
#
#
# def test_quantify_read_barcode_map_files():
#     outMap = "/nfs/home/ansintsova/wishtags/tags.txt"
#
#     barcode_2_pos, barcode_2_abundance = quantify_barcodes.quantify_read_barcode_map_files(outMap)
#     print(barcode_2_pos)
#
# def test_quantify_wishtags():
#     in_fq_file = "/nfs/home/ansintsova/WISHTAGS_PAN/data/20200806_FS10000965_4_BPC29618-1621/MixA/01/MixA_S1_L001_R1_001.fastq.gz"
#     tnSeq = 'GGAGGTTCACAATGTGGGAGGTCA'
#     cnter = quantify_barcodes.quantify_load_fq_barcodes(in_fq_file, tnSeq, bc2tp2=0, bcLen=40, before=False)
#     max_edit_distance = 2
#     outMap = "/nfs/home/ansintsova/WISHTAGS_PAN/data/tags.txt"
#     barcode_2_pos, barcode_2_abundance = quantify_barcodes.quantify_read_barcode_map_files(outMap)
#     testOutFile = "/nfs/home/ansintsova/WISHTAGS_PAN/data/quantifyExtractTest2.out"
#     t = quantify_barcodes.quantify_extract_annotated_correct(barcode_2_pos, barcode_2_abundance, cnter, testOutFile, max_edit_distance)
#     print(t)


#test_quantify_wishtags()
