import collections
import os
from typing import Optional, List, Union, Tuple
from mbarq.core import Barcode, BarSeqData
from pathlib import Path
import pandas as pd
import shutil
import sys
from mbarq.mbarq_logger import get_logger
import numpy as np
from pybedtools import BedTool
import pyranges as pr
from Bio.Seq import Seq


class Mapper(BarSeqData):

    """

    Class to contain and process sequencing data from a mapping BarSeq experiment


    """
    def __init__(self, sequencing_file: str, barcode_structure: str, name: str = '',
                 genome: str = '',
                 output_dir: str = ".",
                 edit_distance: int = 3
                 ) -> None:
        super().__init__(sequencing_file)
        self.barcode_structure = barcode_structure
        self.output_dir = Path(output_dir)
        self.genome = shutil.copy(genome, self.output_dir/Path(genome).name) if Path(genome).is_file() else ''
        self.blastdb = ''
        self.name = name if name else Path(self.seq_file.strip('.gz')).stem
        self.temp_fasta_file = self.output_dir / f"{self.name}.fasta"
        self.temp_blastn_file = self.output_dir / f"{self.name}.blastn"
        self.map_file: Path = self.output_dir / f"{self.name}.map.csv"
        self.blast_columns: str = "qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sstrand"
        self.positions: pd.DataFrame = pd.DataFrame()
        self.edit_distance: int = edit_distance
        self.logger = get_logger('mapping-log', self.output_dir / f"{self.name}_mapping.log")
        self._validate_input()
        self.logger.info("Initializing Barcode Mapper")

    def _validate_input(self) -> None:
        """
        - Check that either genome fasta or blast-db are provided
        - Check that sequencing file exists

        """
        if not self.blastdb and not self.genome:
            self.logger.error("Blast DB or genome file are needed for mapping")
            sys.exit(1)
        elif not Path(self.seq_file).is_file():
            self.logger.error(f'{self.seq_file} could not be found')
            sys.exit(1)
        elif self.genome and str(self.genome).endswith('.gz'):
            self.logger.error(f'Currently cannot accept genome file in compressed format, '
                              f'please provide uncompressed FASTA')

    def extract_barcodes(self, min_host_bases: int = 20) -> None:
        """

         :param min_host_bases:  minimal length of host sequence to consider

        - Adapted from original written by Hans

        """
        self.logger.info('------------------')
        self.logger.info("Extracting Barcodes")
        total_inserts: int = 0
        inserts_with_tp2: int = 0
        inserts_without_tp2: int = 0
        inserts_with_tp2_but_short_bc: int = 0
        inserts_with_tp2_with_good_bc: int = 0
        inserts_with_tp2_with_good_bc_short_host: int = 0
        for total_inserts, r1 in enumerate(self.stream_seq_file()):
            barcode = Barcode(self.barcode_structure)
            if total_inserts % 100000 == 0:
                self.logger.info(f'Processed {total_inserts + 1} reads')
            if barcode.tn_seq in r1.sequence:
                inserts_with_tp2 += 1
                barcode.extract_barcode_host(r1)
                if barcode.bc_seq:
                    inserts_with_tp2_with_good_bc += 1
                    if not barcode.host or len(barcode.host) < min_host_bases:
                        inserts_with_tp2_with_good_bc_short_host += 1
                    else:
                        self.barcodes.append(barcode)
                else:
                    inserts_with_tp2_but_short_bc += 1
            else:
                inserts_without_tp2 += 1
        self.logger.info('Finished Extraction')
        total_inserts += 1
        self.logger.info(f'Processed {total_inserts} reads')
        self.logger.info(f'\tTotal inserts:\t{total_inserts}\t100.0%')
        self.logger.info(f'\tInserts w transposon:\t{inserts_with_tp2}\t'
                         f'{int(100 * ((inserts_with_tp2 * 100.0) / total_inserts)) / 100.0}%')
        self.logger.info(f'\t\tand w barcode:\t{inserts_with_tp2_with_good_bc}\t'
                         f'{int(100 * ((inserts_with_tp2_with_good_bc * 100.0) / total_inserts)) / 100.0}%')
        self.logger.info(f'\t\t\tand good host sequence:\t{len(self.barcodes)}\t'
                         f'{int(100 * ((len(self.barcodes) * 100.0) / total_inserts)) / 100.0}'
                         f'% --> Used for downstream analysis')
        self.logger.info(f'\t\t\tand short host sequence:\t{inserts_with_tp2_with_good_bc_short_host}\t'
                         f'{int(100 * ((inserts_with_tp2_with_good_bc_short_host * 100.0) / total_inserts)) / 100.0}%')
        self.logger.info(f'\t\tand w/o barcode:\t{inserts_with_tp2_but_short_bc}\t'
                         f'{int(100 * ((inserts_with_tp2_but_short_bc * 100.0) / total_inserts)) / 100.0}%')
        self.logger.info(f'\tInserts w/o transposon:\t{inserts_without_tp2}\t'
                         f'{int(100 * ((inserts_without_tp2 * 100.0) / total_inserts)) / 100.0}%')

    def _dereplicate_barcodes(self):
        """

        - Original written by Hans

        """
        self.logger.info('------------------')
        self.logger.info("Dereplicating barcodes")
        barcode_2_sequences = collections.defaultdict(list)
        for barcode in self.barcodes:
            barcode_2_sequences[barcode.bc_seq].append(barcode.host)
        for barcode_seq, sequences in barcode_2_sequences.items():
            for sequence, cnt in collections.Counter(sequences).most_common():
                barcode = Barcode(self.barcode_structure)
                barcode.bc_seq = barcode_seq
                barcode.host = sequence
                barcode.count = cnt
                yield barcode

    def _write_barcodes_to_fasta(self):
        self.logger.info('------------------')
        self.logger.info('Writing dereplicated barcodes to file')
        with open(self.temp_fasta_file, 'w') as handle:
            tot = 1
            for barcode in self._dereplicate_barcodes():
                handle.write(f'>{tot}_bc_{barcode.bc_seq}_cnt_{barcode.count}\n{barcode.host}\n')
                tot += 1

    def _blast_barcode_host(self, blast_threads: int = 4) -> None:
        """
         Map reads against the  blast database
         Creates database from reference genome if it doesn't already exist)
         :return: None
         - Original written by Hans

         """
        db_return_code = 0
        self.logger.info('------------------')
        self.logger.info('Running BLAST on barcode-associated host sequences')
        if not self.blastdb:
            if self.genome.suffix == '.gz':
                self.blastdb = self.genome.with_suffix("")
                command0 = f'gunzip -c {self.genome} | makeblastdb -in - -dbtype nucl -out {self.blastdb} -title {self.blastdb}'
            else:
                self.blastdb = self.genome
                command0 = f'makeblastdb -in {self.genome} -dbtype nucl'
            self.logger.info('No BLAST DB provided, generating BLAST DB from reference genome.')
            self.logger.info(f"Blast command:\t{command0}")
            db_return_code = self.check_call(command0)
            if db_return_code != 0:
                self.logger.error(f"Failed to create blast DB from {self.genome}. "
                                  f"Return code: {db_return_code}")
                sys.exit(1)
        command = f'blastn -task blastn -db {self.blastdb} -out {self.temp_blastn_file} ' \
                  f'-query {self.temp_fasta_file} ' \
                  f'-outfmt "6 {self.blast_columns}" ' \
                  f'-num_threads {blast_threads}'
        self.logger.info(f'blastn command:\t{command}')
        blast_return_code = self.check_call(command)
        if blast_return_code != 0:
            self.logger.error(f"blastn failed. "
                              f"Return code: {db_return_code}")
            sys.exit(1)
        files_to_remove = [self.blastdb.with_suffix(self.blastdb.suffix + i) for i in ['', '.nhr', '.nin', '.nsq']]
        for file in files_to_remove:
            if file.is_file():
                os.remove(file)


    def _find_most_likely_positions(self, filter_below, perc_primary_location=0.75) -> None:
        """
         Takes in blast file, and provides most likely locations for each barcode
         :param: filter_below
         :return: pd.DataFrame
         """

        self.logger.info('------------------')
        self.logger.info('Filtering blast results to find barcode positions')
        chunks = pd.read_table(self.temp_blastn_file, header=None,
                               names="qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sstrand".split(),
                               usecols="qseqid sseqid pident length sstart send evalue bitscore sstrand".split(),
                               chunksize=1000000)
        df = pd.concat([chunk[(chunk.evalue < 0.1) & (chunk.length > 20)] for chunk in chunks])
        barcode = Barcode(self.barcode_structure)
        insertion_col = 'sstart' if barcode.bc_before_tn else 'send'
        df = df.rename({insertion_col: 'insertion_site'}, axis=1)

        # Get a best hit for each qseqID: group by qseqid, find max bitscore
        best_hits = df.groupby('qseqid').agg({'bitscore': ['max']}).reset_index()
        best_hits.columns = ['qseqid', 'bitscore']
        # Get barcode out of qseqid
        best_hits['barcode'] = best_hits['qseqid'].str.split('_', expand=True)[[2]]
        # Get count out of qseqid
        best_hits['cnt'] = best_hits['qseqid'].str.split('_', expand=True)[[4]].astype(int)
        # Only keep barcodes with max bitscores
        query_best_hits = best_hits.merge(df, how='left', on=['qseqid', 'bitscore'])
        query_best_hits['end'] = query_best_hits['insertion_site'].astype(int) + 5
        self.logger.info("Merge similar positions")
        query_best_hits = query_best_hits.sort_values(['barcode', 'sseqid', 'insertion_site'])
        query_best_hits['Group'] = ((query_best_hits.end.rolling(window=2, min_periods=1).min()
                                     - query_best_hits['insertion_site'].rolling(window=2,
                                                                                 min_periods=1).max()) < 0).cumsum()
        query_best_hits['Group'] = query_best_hits.barcode + "_" + query_best_hits.sseqid.astype(
            str) + "_" + query_best_hits.Group.astype(str)
        cnt = query_best_hits.groupby(['Group']).agg({'cnt': ['sum']}).reset_index()
        cnt.columns = ['Group', 'total_count']
        loc = query_best_hits.loc[query_best_hits.groupby(['Group'])['cnt'].idxmax()]
        loc = loc.merge(cnt, on=['Group'])
        total_counts = loc[['barcode', 'sseqid', 'insertion_site', 'sstrand', 'total_count']].copy()
        self.logger.info("Calculate proportion of reads per position")
        total_counts['prop_read_per_position'] = total_counts['total_count'] / total_counts.groupby('barcode')[
            'total_count'].transform('sum')
        likely_positions = total_counts[total_counts['prop_read_per_position'] > perc_primary_location]
        likely_multimappers = (total_counts[(total_counts['prop_read_per_position'] < perc_primary_location)
                                            & (total_counts.total_count > filter_below)]
                               .barcode.nunique())
        self.logger.info(f"Estimated # of multimappers removed: {likely_multimappers}")
        self.logger.info(f"Filtering out barcodes supported by less than {filter_below} reads")
        likely_positions = likely_positions[likely_positions.total_count > filter_below]
        self.positions = likely_positions[
            ['barcode', 'sseqid', 'sstrand', 'insertion_site', 'total_count', 'prop_read_per_position']]

    def _merge_colliding_barcodes(self) -> None:

        """

        Takes data frame of positions, and merges colliding barcodes

        """
        self.logger.info('------------------')
        self.logger.info(f'Merging barcodes with edit distance of less than {self.edit_distance} mapped to '
                         f'the same positions (+/- 5 nt)')
        pps = self.positions[['sseqid', 'insertion_site', 'sstrand', 'barcode', 'total_count']]
        if pps.empty:
            self.logger.error('No mapped barcodes found. Try lowering the filter_below cutoff')
            sys.exit(1)
        positions_sorted = (pps.groupby('sseqid')
                            .apply(pd.DataFrame.sort_values, 'insertion_site')
                            .drop(['sseqid'], axis=1)
                            .reset_index()
                            .drop(['level_1'], axis=1))

        # Get indices for rows with collisions
        collision_index = list(
            positions_sorted[(positions_sorted['insertion_site'].diff() < 5) & (positions_sorted['insertion_site'].diff() >= 0)].index)
        collision_index.extend([i - 1 for i in collision_index if i - 1 not in collision_index])
        collision_index.sort()

        # Barcodes without collisions

        unique = positions_sorted[~positions_sorted.index.isin(collision_index)]

        collisions = positions_sorted.iloc[collision_index]
        if collisions.empty:
            self.positions = unique[['barcode', 'total_count', 'insertion_site', 'sseqid', 'sstrand']]
        else:
            def row_to_barcode(structure, row):
                bc = Barcode(structure)
                bc.bc_seq = row.barcode
                bc.chr = row.sseqid
                bc.start = row.insertion_site
                bc.strand = row.sstrand
                bc.count = row.total_count
                return bc

            bcs = []
            final_bcs = []
            for i, r in collisions.iterrows():
                bcs.append(row_to_barcode(self.barcode_structure, r))

            bc = bcs.pop(0)
            cnt = bc.count
            while len(bcs) > 0:
                bc2 = bcs.pop(0)
                if bc.chr != bc2.chr or abs(bc.start - bc2.start) > 5:
                    bc.count = cnt
                    final_bcs.append(bc)
                    bc = bc2
                    cnt = bc2.count
                else:
                    if bc.editdistance(bc2) > self.edit_distance:
                        bc.count = cnt
                        final_bcs.append(bc)
                        bc = bc2
                        cnt = bc2.count
                    else:
                        cnt += bc2.count
                        if bc.count < bc2.count:
                            bc = bc2
            if bc not in final_bcs:
                final_bcs.append(bc)
            resolved_collisions = pd.DataFrame([[bc.chr, bc.start, bc.strand,
                                                 bc.bc_seq, bc.count] for bc in final_bcs],
                                               columns=['sseqid', 'insertion_site', 'sstrand', 'barcode', 'total_count'])
            self.positions = pd.concat([unique, resolved_collisions])[['barcode', 'total_count', 'insertion_site',
                                                                       'sseqid', 'sstrand']]
        self.logger.info('------------------')
        self.logger.info("Finished mapping barcodes")
        self.logger.info(f"Final number of barcodes found: {self.positions.barcode.nunique()}")
        self.logger.info(f"Mean number of reads per barcode: {self.positions.total_count.mean()}")

    def map_insertions(self, min_host_bases: int = 20, filter_below: int = 100, blast_threads: int =4,
                       no_blast: bool = False):

        """
        :param min_host_bases: Minimal length of host sequence to consider
        :param filter_below: Read count threshold for keeping barcodes
        :param no_blast: Whether to run the blast command or not
        0. If no_blast == True, look for existing blast output file.
        1. Else Extract barcodes from the sequencing file using `extract_barcodes`, write barcode to fasta, blast each host sequence agains the genome
        2. Find most likely positions
        3. Merge colliding barcodes
        4. Write out final positions to csv

        """

        if no_blast:
            if not self.temp_blastn_file.is_file():
                self.logger.error("--no_blast option specified, but no BLAST output file found")
                sys.exit(1)
        else:
            self.extract_barcodes(min_host_bases=min_host_bases)
            self._write_barcodes_to_fasta()
            self._blast_barcode_host(blast_threads)
        self._find_most_likely_positions(filter_below=filter_below)
        self._merge_colliding_barcodes()
        rename_dict = {'sseqid': 'chr',
                       'sstrand': 'strand',
                       'total_count': 'abundance_in_mapping_library'}
        self.positions = self.positions.rename(rename_dict, axis=1)
        self.logger.info('Writing transposon insertion sites to file')
        self.positions.to_csv(self.map_file, index=False)
        self.logger.info("Removing intermediate files")
        os.remove(self.temp_fasta_file)
        os.remove(self.temp_blastn_file)
        self.logger.info("Done")


class AnnotatedMap:
    def __init__(self, map_file: str,
                 annotation_file: str,
                 feature_type: str,
                 identifiers: Tuple[str, ...],
                 name: str = '',
                 output_dir: str = ".",
                 positions: pd.DataFrame = pd.DataFrame()
                 ) -> None:
        self.map_file: Union[str, Path] = map_file
        if name:
            self.name: str = name
        else:
            self.name = Path(self.map_file).stem
        self.output_dir: Path = Path(output_dir)
        self.logger = get_logger('annotate_map-log', self.output_dir / f"{self.name}_annotate_map.log")
        self.annotations: str = annotation_file
        self.feature_type: str = feature_type
        self.identifiers: List[str, ...] = list(identifiers)
        self.annotated_map_file: Path = self.output_dir / f'{self.name}.annotated.csv'
        self.positions: pd.DataFrame = positions if not positions.empty else pd.read_csv(self.map_file)
        self.temp_bed_file: Path = self.output_dir / f"{self.name}.bed"

    def _find_annotation_overlaps(self, intersect=True):

        """

        """

        def min_max(isite, start, end):
            return round((isite - start) / (end - start), 2)

        # self.logger.info('------------------')
        # self.logger.info('Annotating mapped barcodes')
        barcode_cols = ['chr', 'insertion_site', 'barcode', 'abundance_in_mapping_library']
        bed_map = self.positions[barcode_cols].copy()
        bed_map['bc_start'] = bed_map['insertion_site'] - 1
        barcode_cols = ['chr', 'bc_start', 'insertion_site', 'barcode', 'abundance_in_mapping_library']
        bed_map = bed_map[barcode_cols]
        barcodes_bedtool = BedTool.from_dataframe(bed_map).sort()
        gff = pr.read_gff3(self.annotations).as_df()
        gff['Chromosome'] = gff['Chromosome'].astype(str)
        gff['Start'] = gff['Start'] + 1
        gff_cols = ['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame']
        gff_short = gff[gff.Feature == self.feature_type][gff_cols]
        genes_bedtool = BedTool.from_dataframe(gff_short).sort()
        names = barcode_cols + gff_cols + ['distance_to_feature']
        barcode_annotation = barcodes_bedtool.closest(genes_bedtool, D='b').to_dataframe(names=names)
        if intersect == True:
            barcode_annotation = barcode_annotation[barcode_annotation['distance_to_feature'] == 0]
        # calculate percentile
        barcode_annotation['Chromosome'] = barcode_annotation['Chromosome'].astype(str)
        barcode_annotation = barcode_annotation.merge(gff[gff_cols + self.identifiers], on=gff_cols,
                                                      how='left')
        barcode_annotation['percentile'] = np.where(barcode_annotation.Strand.isin(['+', 'plus']),
                                                    min_max(barcode_annotation.insertion_site, barcode_annotation.Start,
                                                            barcode_annotation.End),
                                                    min_max(barcode_annotation.insertion_site, barcode_annotation.End,
                                                            barcode_annotation.Start))
        barcode_annotation['percentile'] = barcode_annotation['percentile'].apply(
            lambda x: np.nan if x <= 0.0 or x >= 1.0 else x)
        final_columns = ['barcode', 'chr', 'insertion_site', 'abundance_in_mapping_library',
                         'Start', 'End', 'Strand'] + self.identifiers + ['distance_to_feature', 'percentile']
        self.annotated_positions = (barcode_annotation[final_columns]
                                    .rename({'Start': 'gene_start', 'End': 'gene_end',
                                             'Strand': 'gene_strand'}, axis=1))
        self.annotated_positions.to_csv(self.annotated_map_file, index=False)

    def _validate_annotations(self) -> None:

        """
        Simple validation on the gff file
        only if it's a non-empty file with feature type

        """
        if not Path(self.annotations).is_file():
            self.logger.error(f'Annotation file "{self.annotations}" not found. Check the path.')
            self.logger.error(f'To rerun annotation on {self.map_file} run ".... "')
            sys.exit(1)
        else:
            feature_types = pd.read_table(self.annotations, comment='#', usecols=[2]).iloc[:, 0].unique()
            if self.feature_type not in feature_types:
                # todo assumes gff format, might want to generalize later
                self.logger.error(f'Feature type "{self.feature_type}" not found in {self.annotations}.')
                self.logger.error(f'Feature types available: {", ".join(feature_types)}')
                sys.exit(1)

    def add_rev_complement(self):
        self.logger.info(f"Adding barcode reverse complements to {self.annotated_map_file}")
        self.annotated_positions['revcomp_barcode'] = self.annotated_positions.barcode.apply(
            lambda x: str(Seq(x).reverse_complement()))
        self.logger.info(f"writing to {self.annotated_map_file}")
        self.annotated_positions.to_csv(self.annotated_map_file, index=False)

    def annotate(self, intersect=True):
        if self.positions.empty:
            self.logger.error(f'Found no positions to annotate. Is {self.map_file} empty?')
            sys.exit(1)
        self._validate_annotations()
        self._find_annotation_overlaps(intersect)
        self.logger.info("Done!")



