import collections
import os
from typing import Optional, List, Counter
from mbarq.core import Barcode, BarSeqData
from pathlib import Path
import pandas as pd
import sys
from mbarq.mbarq_logger import get_logger


class BarcodeCounter(BarSeqData):

    def __init__(self, sequencing_file: str,
                 barcode_structure: str,
                 name: str = '',
                 mapping_file: str = '',
                 output_dir: str = ".",
                 edit_distance: int = 3
                 ) -> None:
        super().__init__(sequencing_file)
        self.barcode_structure = barcode_structure
        self.output_dir: Path = Path(output_dir)
        self.name = name if name else Path(self.seq_file.strip('gz')).stem
        self.map_file: str = mapping_file
        self.merged: bool = False
        self.barcode_counter: Counter = collections.Counter()
        self.edit_distance: int = edit_distance
        self.annotated_cnts: pd.DataFrame = pd.DataFrame()
        self.counts_file: Path = self.output_dir/ f"{self.name}_mbarq_counts.csv"
        self.logger = get_logger('counting-log', self.output_dir / f"{self.name}_counting.log")
        self.validate_input()
        self.logger.info("Initializing Counter")

    def _extract_barcodes(self) -> None:
        """
        Count barcodes in a sequencing file (fasta/fastq)
        """
        self.logger.info('------------------')
        self.logger.info("Extracting Barcodes")
        self.logger.info("Counting reads with and without transposon sequence")
        total_inserts: int = 0
        with_tp2: int = 0
        without_tp2: int = 0
        with_tp_but_short: int = 0
        for total_inserts, r1 in enumerate(self.stream_seq_file(), 1):
            barcode = Barcode(self.barcode_structure)
            if total_inserts % 1000000 == 0:
                self.logger.info(f'\tReads processed:\t{total_inserts}')
            if barcode.tn_seq in r1.sequence:
                barcode.extract_barcode_host(r1)
                if not barcode.bc_seq:
                    with_tp_but_short += 1
                else:
                    with_tp2 += 1
                    self.barcode_counter[barcode.bc_seq] += 1
            else:
                without_tp2 += 1
        self.logger.info(f'Reads processed:\t{total_inserts}')
        self.logger.info(f'FastA/Q Stats:')
        self.logger.info(f'\tTotal Reads:\t{total_inserts}')
        self.logger.info(f'\tReads w transposon and good barcode:\t{with_tp2}')
        self.logger.info(f'\tReads w transposon but short barcode:\t{with_tp_but_short}')
        self.logger.info(f'\tReads w/o transposon:\t{without_tp2}')
        self.logger.info(f'\t% Reads w/o transposon: \t {without_tp2/total_inserts*100}')

    def _merge_similar(self):

        for sequence, count in self.barcode_counter.most_common():
            b = Barcode(sequence=sequence)
            b.count = count
            self.barcodes.append(b)
        for rare_barcode in self.barcodes[::-1]:
            idx = self.barcodes.index(rare_barcode)
            for common_barcode in self.barcodes[:idx]:
                if common_barcode.editdistance(rare_barcode) < self.edit_distance:
                    common_barcode.count += rare_barcode.count
                    self.barcodes.remove(rare_barcode)
                    break
        self.barcode_counter = collections.Counter({barcode.bc_seq: barcode.count for barcode in self.barcodes})
        self.merged = True

    def _annotate_barcodes(self) -> None:
        assert self.barcode_counter is not None
        self.logger.info(f'Barcodes with edit distance of less than {self.edit_distance} '
                         f'have been merged: {self.merged}')
        cnts_df = (pd.DataFrame.from_dict(self.barcode_counter, orient='index')
                   .reset_index())
        cnts_df.columns = ['barcode', 'barcode_count']
        cnts_df = cnts_df[cnts_df['barcode_count'] > 1]

        if cnts_df.empty:
            self.logger.error('No barcodes with counts > 1 found')
            sys.exit(1)
        self.logger.info(f'Number of unique barcodes to annotate: {cnts_df.barcode.nunique()}')

        if not self.map_file:
            self.logger.info('No mapping file found, skipping the annotation')
            self.annotated_cnts = cnts_df
        else:
            self.logger.info(f'Mapping file used: {self.map_file}')
            bc_annotations = pd.read_csv(self.map_file)
            to_keep = bc_annotations.columns
            if 'barcode' not in bc_annotations.columns:
                self.logger.error('No column "barcode" found in the mapping/annotation file')
                sys.exit(1)
            if 'barcode_count' in bc_annotations.columns:
                self.logger.error('"barcode_count" column already present in the '
                                  'mapping/annoations file. Please rename to avoid confusion')
                sys.exit(1)
            self.logger.info(f'Columns found in the mapping/annotations file: '
                             f'{", ".join(bc_annotations.columns)}')
            self.annotated_cnts = cnts_df.merge(bc_annotations, how='left', on='barcode')
            filter_col = to_keep[1]
            with_ids = self.annotated_cnts[self.annotated_cnts[filter_col].notnull()]
            self.logger.info(f'Number of annotated barcodes: {with_ids.shape[0]}')
            no_ids = self.annotated_cnts[self.annotated_cnts[filter_col].isna()]
            self.logger.info(f'Number of unannotated barcodes:{no_ids.shape[0]}')

    def _write_counts_file(self) -> None:
        self.annotated_cnts.to_csv(self.counts_file)

    def count_barcodes(self):
        self.logger.info("1")
        self._extract_barcodes()
        self.logger.info("2")
        self._merge_similar()
        self.logger.info("3")
        self._annotate_barcodes()
        self.logger.info(f'Writting final counts to {self.counts_file}')
        self._write_counts_file()
