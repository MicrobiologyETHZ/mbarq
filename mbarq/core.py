from typing import Optional, Tuple, List
from pathlib import Path
import gzip
import Bio.SeqIO.QualityIO
from Bio.Seq import Seq
import Bio.SeqIO.FastaIO as FastaIO
from typing import Generator
import logging
import subprocess
import re


class InputFileError(Exception):
    pass


class FastA:
    '''
    Standard data container for fasta sequences
    '''
    __slots__ = ['header', 'sequence']

    def __init__(self, header: str, sequence: str) -> None:
        self.header = header
        self.sequence = sequence


class Barcode:
    def __init__(self, structure='', sequence='', host=''):
        self.structure: str = structure
        self.bc_seq: str = sequence
        self.host: str = host
        self.bc_len: int
        self.tn_seq: str
        self.count: int = -1
        self.bc_before_tn: bool
        self.len_spacer: int
        # In theory these are optional
        self.start: Optional[int] = None  # todo don't need these, need insertion site
        self.end: Optional[int] = None
        self.chr: Optional[str] = None
        self.strand: Optional[str] = None
        self.multimap: Optional[bool] = None
        self.identifiers: Optional[List[str]] = None
        if self.structure:
            self._parse_structure()
        if not self.structure and not self.bc_seq:
            raise ValueError("Please provide either structure or sequence")

    # def _parse_structure_old(self):
    #     self.tn_seq = self.structure.split(':')[0]
    #     self.len_spacer = int(self.structure.split(':')[2])
    #     self.bc_len = int(self.structure.split(':')[1])
    #     self.bc_before_tn = True if self.structure.split(':')[3] == 'before' else False

    def _parse_structure(self):
        try:
            self.tn_seq = re.findall('[ACGT]+', self.structure)[0]
            bc_len = re.findall('B(\\d+)', self.structure)[0]
            spacer = re.findall('N(\\d+)', self.structure)
            self.len_spacer = int(spacer[0]) if spacer else 0
            self.bc_before_tn = self.structure.index(self.tn_seq) > self.structure.index(bc_len)
            self.bc_len = int(bc_len)


        except IndexError:
            raise 'Could not transposon structure provided'



    def editdistance(self, other_barcode: "Barcode") -> int:
        '''
        Calculate the edit distance between 2 sequences with identical length.
        Will throw an error if the length of both sequences differs

        :return:
        '''
        if not self.bc_seq or not other_barcode.bc_seq:
            raise Exception("Could not find sequence for one or more barcodes")
        if self.bc_seq == other_barcode.bc_seq:
            return 0
        if len(self.bc_seq) != len(other_barcode.bc_seq):
            raise Exception(
                f'{self.bc_seq} and {other_barcode.bc_seq} have different length. Edit distance can be computed on same length sequences only.')
        dist = 0
        for letter1, letter2 in zip(self.bc_seq, other_barcode.bc_seq):
            if letter1 != letter2:
                dist += 1
        return dist

    def extract_barcode_host(self, r1: FastA) -> None:
        '''
       Extract barcode and host sequence from read with tp2.
       Return (None, None) if the barcode sequence is not complete (17bp)
        :param r1:
        :return:

        -----|BARCODE|----------|TN end sequence (tp2)|---Host------
        -----|-bc_len-|--len_spacer--|---------tn_seq---------|-------------
        -----|-17bp--|---13bp---|---------15bp--------|----?--------
         ---(-30)---(-13)-------(0)---------------------------------
        '''
        splits: List[str] = r1.sequence.split(self.tn_seq)  # check that tn in sequence?
        self.bc_seq = ''
        self.host = ''
        if self.bc_before_tn:
            bc_start = -(self.bc_len + self.len_spacer)
            bc_end = None if self.len_spacer == 0 else -self.len_spacer
            bc_seq = splits[0]
            host_seq = splits[1]
        else:
            bc_start = self.len_spacer
            bc_end = self.len_spacer + self.bc_len
            bc_seq = splits[1]
            host_seq = splits[0]

        if len(bc_seq) >= self.len_spacer + self.bc_len:
            self.bc_seq = bc_seq[bc_start:bc_end]
            self.host = host_seq

    def __repr__(self):
        if self.bc_seq:
            return f"{self.bc_seq}: {self.count}"
        else:
            return f"Barcode({self.structure})"


class BarSeqData:
    def __init__(self, sequencing_file: str, annotation_file: str = '', ) -> None:
        self.seq_file = str(sequencing_file)
        self.annotations = Path(annotation_file)
        self.barcodes: List[Barcode] = []

    def _validate_input(self) -> None:
        if not Path(self.seq_file).is_file():
            raise InputFileError(f'{self.seq_file} could not be found')

    def stream_fastq_file(self):
        self._validate_input()
        if self.seq_file.endswith('fq.gz') or self.seq_file.endswith('fastq.gz'):
            with gzip.open(self.seq_file, 'rt') as handle:
                for record in Bio.SeqIO.parse(handle, 'fastq'):
                    yield record
        elif self.seq_file.endswith('fq') or self.seq_file.endswith('fastq'):
            with open(self.seq_file) as handle:
                for record in Bio.SeqIO.parse(handle, 'fastq'):
                    yield record
        else:
            raise Exception(f'{self.seq_file} not a FASTQ file.')

    def stream_seq_file(self) -> Generator[FastA, None, None]:
        """
        Read a fastq file either gzipped or not and return it as a stream of tuples
        (Header, Sequence, Quality)
        :param infile:
        :return: Generator[FastA, None, None]
        """
        self._validate_input()
        if self.seq_file.endswith('fq.gz') or self.seq_file.endswith('fastq.gz'):
            with gzip.open(self.seq_file, 'rt') as handle:
                for header, sequence, qual in Bio.SeqIO.QualityIO.FastqGeneralIterator(handle):
                    yield FastA(header, sequence)
        elif self.seq_file.endswith('fq') or self.seq_file.endswith('fastq'):
            with open(self.seq_file) as handle:
                for header, sequence, qual in Bio.SeqIO.QualityIO.FastqGeneralIterator(handle):
                    yield FastA(header, sequence)
        elif self.seq_file.endswith('fasta.gz') or self.seq_file.endswith('fa.gz'):
            with gzip.open(self.seq_file, 'rt') as handle:
                for (header, sequence) in FastaIO.SimpleFastaParser(handle):
                    yield FastA(header, sequence)
        elif self.seq_file.endswith('fasta') or self.seq_file.endswith('fa'):
            with open(self.seq_file) as handle:
                for (header, sequence) in FastaIO.SimpleFastaParser(handle):
                    yield FastA(header, sequence)
        else:
            raise Exception(f'{self.seq_file} not a sequence file.')

    def check_call(self, command) -> int:
        """
        Simple wrapper to execute check_call and catch exceptions
        :param command:
        :return:
        """
        returncode = -1
        try:
            returncode: int = subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as e:
            raise Exception(e)
        #if returncode != 0:
            #raise (Exception(f'Command: {command} failed with return code {returncode}'))
        return returncode

    def __repr__(self) -> str:
        return f"{self.seq_file}"


