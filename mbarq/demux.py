from mbarq.core import BarSeqData
from Bio.Seq import Seq
from pathlib import Path
from Bio import SeqIO
import sys
from typing import Union
import subprocess

def demux_barseq(sequencing_file: Union[str, Path],
                demux_bc_file: Union[str, Path], out_dir: Union[str, Path] ='.',
                name: str = '',  rc: bool = False, sep='\t'):
    """
    Takes in multiplexed fastq file and demux_bc_file
    Example demux_bc_file (tsv): H1 ACCGTG
                           H2 AGGCCT

    """
    index_2_sample = {}
    len_barcodes = set()
    with open(demux_bc_file, 'r') as fh:
        for line in fh.readlines():
            line = line.strip()
            if line:
                samplename, index_barcode = line.split(sep)
                if name:
                    samplename = name + "_" + samplename
                if rc:
                    index_barcode = str(Seq(index_barcode).reverse_complement())
                index_2_sample[index_barcode] = samplename
                len_barcodes.add(len(index_barcode))
    if len(len_barcodes) > 1:
        print("Not all barcodes the same length")
        sys.exit(1)
    else:
        len_index = list(len_barcodes)[0]
    index_2_writer = {}
    for index_barcode, samplename in index_2_sample.items():
        index_2_writer[index_barcode] = open(Path(out_dir) / f"{samplename}.fastq", 'w')
    good_index = 0
    barseq = BarSeqData(sequencing_file=sequencing_file)
    for total_inserts, r1 in enumerate(barseq.stream_fastq_file(), start=1):
        index = r1.seq[0:len_index]
        if index in index_2_sample:
            good_index += 1
            SeqIO.write(r1, index_2_writer[index], "fastq")
        if total_inserts % 100000 == 0:
            print(f'Good/Total\t{good_index}/{total_inserts}')
    for index, writer in index_2_writer.items():
        writer.close()
    for samplename in index_2_sample.values():
        print(f'Compressing {samplename}.fastq')
        subprocess.check_call(['gzip', Path(out_dir)/f"{samplename}.fastq"])


if __name__ == "__main__":
    root = Path("/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/nguyenb/tnseq/scratch/deutschbauer/fastq")
    sequencing_file = root/"SB2B_ML5_set1_sub01.fastq.gz"
    barcode_file = root/"demux_barcodes.txt"
    out_dir = root/"demux"
    demux_barseq(sequencing_file, barcode_file, out_dir)