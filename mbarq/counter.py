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
                 edit_distance: int = 2,
                 rev_complement: bool = False,
                 ) -> None:
        super().__init__(sequencing_file)
        self.output_dir: Path = Path(output_dir)
        self.name = name if name else Path(self.seq_file.strip('.gz')).stem
        self.logger = get_logger('counting-log', self.output_dir / f"{self.name}_counting.log")
        self.logger.info("Initializing Barcode Counter")
        self.barcode_structure = barcode_structure
        self.map_file: str = mapping_file
        self.merged: bool = False
        self.barcode_counter: Counter = collections.Counter()
        self.edit_distance: int = edit_distance
        self.rev_complement: bool = rev_complement
        self.annotated_cnts: pd.DataFrame = pd.DataFrame()
        # Important to keep _mbarq in the final count file name -> used to get sample IDs during merge
        self.counts_file: Path = self.output_dir / f"{self.name}_mbarq_counts.csv"
        self.bc_annotations: pd.DataFrame = self._validate_annotations()
        self._validate_input()

    def _validate_annotations(self) -> pd.DataFrame:
        """

        Simple validation on the mapping file
        :return: pd.DataFrame

        """
        if not self.map_file:
            self.logger.info(f'No barcode mapping/annotations file provided. Annotation step will be skipped')
            return pd.DataFrame(columns=['barcode'])
        else:
            self.logger.info(f'Annotation file provided: {self.map_file}')
            bc_annotations = pd.read_csv(self.map_file)
            if 'barcode_count' in bc_annotations.columns:
                self.logger.error('"barcode_count" column already present in the '
                                  'mapping/annotation file. Please rename to avoid confusion')
                sys.exit(1)
            if self.rev_complement:
                if 'revcomp_barcode' not in bc_annotations.columns:
                    self.logger.error('No column "revcomp_barcode" found in the mapping/annotation file')
                    sys.exit(1)
                elif 'barcode' in bc_annotations.columns:
                    bc_annotations = bc_annotations.drop(['barcode'], axis=1)
                return bc_annotations.rename({'revcomp_barcode': 'barcode'}, axis=1)
            if 'barcode' not in bc_annotations.columns:
                self.logger.error('No column "barcode" found in the mapping/annotation file')
                sys.exit(1)
            return bc_annotations

    def _extract_barcodes(self) -> None:
        """
        Count barcodes in a sequencing file (fasta/fastq)

        """
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

    def _merge_similar(self) -> None:
        """

        If edit_distance > 0, merge barcodes with edit distance < or equal to edit_distance
        :return: None

        """
        for sequence, count in self.barcode_counter.most_common():
            b = Barcode(sequence=sequence)
            b.count = count
            self.barcodes.append(b)

        if self.edit_distance < 1:
            self.logger.info(f'Not merging similar barcodes')
            self.logger.info(f'Number of unique barcodes: {len(self.barcodes)}')
            self.merged = False
        else:
            self.logger.info(f'Merging barcodes with edit distance <= {self.edit_distance}')
            self.logger.info(f"Number of barcodes to process: {len(self.barcodes)}")
            processed_barcodes = 0
            for rare_barcode in self.barcodes[::-1]:
                idx = self.barcodes.index(rare_barcode)
                processed_barcodes += 1
                if processed_barcodes % 1000 == 0:
                    self.logger.info(f'\tBarcodes processed:\t{processed_barcodes}')
                for common_barcode in self.barcodes[:idx]:
                    if common_barcode.editdistance(rare_barcode) <= self.edit_distance:
                        if rare_barcode.bc_seq in self.bc_annotations['barcode'].values and not common_barcode.bc_seq in self.bc_annotations['barcode'].values:
                            common_barcode.bc_seq = rare_barcode.bc_seq
                        common_barcode.count += rare_barcode.count
                        self.barcodes.remove(rare_barcode)
                        break
            self.merged = True
        self.barcode_counter = collections.Counter({barcode.bc_seq: barcode.count for barcode in self.barcodes})

    def _annotate_barcodes(self, filter_low=False, annotated_only=False) -> None:
        """
        If annotation is available, annotate the counted barcodes
        :return: None
        """
        assert self.barcode_counter is not None
        self.logger.info(f'Barcodes with edit distance of less than {self.edit_distance} '
                         f'have been merged: {self.merged}')
        cnts_df = (pd.DataFrame.from_dict(self.barcode_counter, orient='index')
                   .reset_index())
        # Important to keep these column names, used for merging
        cnts_df.columns = ['barcode', 'barcode_count']
        if filter_low:
            cnts_df = cnts_df[cnts_df['barcode_count'] > int(filter_low)]
        if cnts_df.empty:
            self.logger.error('No barcodes with counts > 0 found')
            sys.exit(1)
        self.logger.info(f'Number of unique barcodes to annotate: {cnts_df.barcode.nunique()}')
        if self.bc_annotations.empty:
            self.logger.info('No mapping file provided, skipping the annotation')
            self.annotated_cnts = cnts_df
        else:
            self.logger.info(f'Annotating barcodes')
            self.logger.info(f'Columns found in the mapping/annotations file: '
                             f'{", ".join(self.bc_annotations.columns)}')
            merge_strategy = 'inner' if annotated_only else 'left'
            self.annotated_cnts = cnts_df.merge(self.bc_annotations, how=merge_strategy, on='barcode').drop_duplicates()
            filter_col = self.bc_annotations.columns[1]
            with_ids = self.annotated_cnts[self.annotated_cnts[filter_col].notnull()]
            self.logger.info(f'Number of annotated barcodes: {with_ids.shape[0]}')
            no_ids = self.annotated_cnts[self.annotated_cnts[filter_col].isna()]
            self.logger.info(f'Number of unannotated barcodes:{no_ids.shape[0]}')

    def _write_counts_file(self) -> None:
        self.annotated_cnts.to_csv(self.counts_file, index=False)

    def count_barcodes(self, filter_low=False, annotated_only=False):
        self.logger.info('------------------')
        self.logger.info("Step 1: Counting ")
        self.logger.info('------------------')
        self._extract_barcodes()
        self.logger.info('------------------')
        self.logger.info("Step 2: Merging")
        self.logger.info('------------------')
        self._merge_similar()
        self.logger.info('------------------')
        self.logger.info("Step 3: Annotating")
        self.logger.info('------------------')
        self._annotate_barcodes(filter_low, annotated_only)
        self.logger.info("Finished!")
        self.logger.info(f'Writing final counts to {self.counts_file}')
        self._write_counts_file()
