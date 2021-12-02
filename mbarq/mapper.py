import logging
import collections
import os
from typing import Optional, List
from mbarq.core import Barcode, BarSeqData
from pathlib import Path
import pandas as pd
import sys
from mbarq.mbarq_logger import get_logger


class Mapper(BarSeqData):

    def __init__(self, sequencing_file: str, barcode_structure: str, name: str = '',
                 genome: str = '', db: str = '',
                 annotation_file: str = '',
                 identifiers: Optional[List[str]] = None,
                 output_dir: str = ".",
                 edit_distance: int = 3
                 ) -> None:
        super().__init__(sequencing_file, annotation_file)
        self.identifiers = identifiers
        self.barcode_structure = barcode_structure
        self.genome = genome
        self.blastdb = db
        self.output_dir = Path(output_dir)
        self.name = name if name else Path(self.seq_file.strip('gz')).stem
        self.temp_fasta_file = self.output_dir / f"{self.name}.fasta"
        self.temp_bed_file = self.output_dir / f"{self.name}.bed"
        self.temp_bed_results_file = self.output_dir / f"{self.name}.bed.intersect.tab"
        self.temp_bed_closest_file = self.output_dir / f"{self.name}.bed.closest.tab"
        self.temp_blastn_file = self.output_dir / f"{self.name}.blastn"
        self.map_file: Path = self.output_dir / f"{self.name}.map.csv"
        self.blast_columns: str = "qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sstrand"
        self.positions: pd.DataFrame = pd.DataFrame()
        self.annotated_positions: pd.DataFrame = pd.DataFrame()
        self.edit_distance: int = edit_distance
        self.logger = get_logger('mapping-log', self.output_dir / f"{self.name}_mapping.log")
        self.validate_input()
        self.logger.info("Initializing Mapper")

    def validate_input(self) -> None:
        if not self.blastdb and not self.genome:
            self.logger.error("Blast DB or genome file are needed for mapping")
            sys.exit(1)
        elif not Path(self.seq_file).is_file():
            self.logger.error(f'{self.seq_file} could not be found')
            sys.exit(1)

    def extract_barcodes(self, min_host_bases: int = 20) -> None:
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
                self.logger.info(f'Processed {total_inserts} reads')
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

    def _blast_barcode_host(self, blast_threads: int = 4):
        """
         Map reads against the  blast database (creates database from reference genome
          if doesn't already exist)
         :param temp_fasta_file:
         :param temp_blastn_file:
         :param genome:
         :param blastdb:
         :param blast_threads:
         :param logger:
         :return: None

         """
        db_return_code = 0
        self.logger.info('------------------')
        self.logger.info('Running BLAST on barcode-associated host sequences')
        if not self.blastdb:
            command0 = f'makeblastdb -in {self.genome} -dbtype nucl'
            self.logger.info('No BLAST DB provided, generating BLAST DB from reference genome.')
            self.logger.info(f"Blast command:\t{command0}")
            db_return_code = self.check_call(command0)
            if db_return_code != 0:
                self.logger.error(f"Failed to create blast DB from {self.genome}. "
                                  f"Return code: {db_return_code}")
                sys.exit(1)
            self.blastdb = self.genome
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

    def _find_most_likely_positions(self, filter_below) -> None:
        """
         Takes in blast file, and provides most likely locations for each barcode
         :param: blast_file
         :param: filter_below
         :param: logger
         :return: pd.DataFrame
         """
        self.logger.info('------------------')
        self.logger.info('Filtering blast results to find barcode positions')
        df = pd.read_table(self.temp_blastn_file, header=None)
        df.columns = self.blast_columns.split()
        # Filter out spurious hits
        df = df[(df.evalue < 0.001) & (df.pident > 95) & (df.length > 40)]
        # Get a best hit for each qseqID: group by qseqid, find max bitscore
        best_hits = df.groupby('qseqid').agg({'bitscore': ['max']}).reset_index()
        best_hits.columns = ['qseqid', 'bitscore']
        # Get barcode out of qseqid
        best_hits['barcode'] = best_hits['qseqid'].str.split('_', expand=True)[[2]]
        # Get count out of qseqid
        best_hits['cnt'] = best_hits['qseqid'].str.split('_', expand=True)[[4]].astype(int)
        # Note: Total counts are calculated with cnt 1 included,
        # but low counts are filtered out right after
        total_count = best_hits.groupby('barcode').cnt.sum().reset_index()
        total_count.columns = ['barcode', 'total_count']
        best_hits = best_hits.merge(total_count, how='left', on='barcode')
        self.logger.info(f"Filtering out barcodes supported by less than {filter_below} reads")
        best_hits = best_hits[best_hits.cnt > filter_below]
        self.logger.info(f"Number of barcodes found: {best_hits.barcode.nunique()}")
        if best_hits.barcode.nunique() == 0:
            self.logger.error('No mapped barcodes found. Try lowering the filter_below cutoff')
            sys.exit(1)
        # Create best hits data frame by merging best_hits with other columns from blast file
        # There still could be multiple hits for each qseqid, if they have the same blast score
        query_best_hits = best_hits.merge(df, how='left', on=['qseqid', 'bitscore'])
        # Check for multimapping
        multimap = (query_best_hits.groupby(['barcode']).sstart.std(ddof=0) > 5).reset_index().rename(
            {'sstart': 'multimap'},
            axis=1)
        query_best_hits = query_best_hits.merge(multimap, on='barcode')
        # For each barcode select the position supported by most reads
        query_best_hits = query_best_hits.sort_values(['barcode', 'cnt'], ascending=False)
        query_best_hits['rank'] = query_best_hits.groupby(['barcode']).cumcount()
        query_best_hits = query_best_hits[query_best_hits['rank'] == 0].copy()
        query_best_hits.drop('rank', axis=1, inplace=True)
        self.positions = query_best_hits


    def _merge_colliding_barcodes(self) -> None:

        """
        Takes data frame of positions, and merges colliding barcodes
        """
        self.logger.info('------------------')
        self.logger.info(f'Merging barcodes with edit distance of less than {self.edit_distance} mapped to '
                         f'the same positions (+/- 5 nt)')
        pps = self.positions[['sseqid', 'sstart', 'sstrand', 'barcode', 'total_count', 'multimap']]
        if pps.empty:
            self.logger.error('No mapped barcodes found. Try lowering the filter_below cutoff')
            sys.exit(1)
        positions_sorted = (pps.groupby('sseqid')
                            .apply(pd.DataFrame.sort_values, 'sstart')
                            .drop(['sseqid'], axis=1)
                            .reset_index()
                            .drop(['level_1'], axis=1))

        # Get indices for rows with collisions
        collision_index = list(
            positions_sorted[(positions_sorted.sstart.diff() < 5) & (positions_sorted.sstart.diff() >= 0)].index)
        collision_index.extend([i - 1 for i in collision_index if i - 1 not in collision_index])
        collision_index.sort()

        # Barcodes without collisions

        unique = positions_sorted[~positions_sorted.index.isin(collision_index)]

        collisions = positions_sorted.iloc[collision_index]

        def row_to_barcode(structure, row):
            bc = Barcode(structure)
            bc.bc_seq = row.barcode
            bc.chr = row.sseqid
            bc.start = row.sstart
            bc.strand = row.sstrand
            bc.multimap = row.multimap
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
                if bc.editdistance(bc2) > 3:
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
                                             bc.bc_seq, bc.count, bc.multimap] for bc in final_bcs],
                                           columns=['sseqid', 'sstart', 'sstrand', 'barcode', 'total_count',
                                                    'multimap'])
        self.positions = pd.concat([unique, resolved_collisions])[['barcode', 'total_count', 'sstart',
                                                                   'sseqid', 'sstrand', 'multimap']]
        self.logger.info('------------------')
        self.logger.info("Finished mapping barcodes")
        self.logger.info(f"Final number of barcodes found: {self.positions.barcode.nunique()}")
        self.logger.info(f"Mean number of reads per barcode: {self.positions.total_count.mean()}")
        self.logger.info(f"Number of barcodes mapping to multiple locations: {self.positions.multimap.sum()}")

    def _find_annotation_overlaps(self):

        """
        Takes output of merge colliding bcs, turns it into bed file, then finds intersections with
        annotation file.
        Generates tab file with the following columns: chr | sstart | gff-info-field | barcode
        :param barcode_map:
        :param bed_file: file path to create tmp bed file needed by bedtools intersect
        :param gff_file:
        :param output_map: tab file with the following columns: chr | sstart | gff-info-field | barcode
        :return: bedtoosl intersect return code
        """
        self.logger.info('------------------')
        self.logger.info('Annotating mapped barcodes. Only annotating barcodes inside/overlapping features of interest')
        assert self.annotations is not None
        bed_map = self.positions.copy().reset_index()
        bed_map['startOffBy1'] = bed_map['sstart'] - 1
        bed_map[['sseqid', 'startOffBy1', 'sstart', 'barcode']].to_csv(self.temp_bed_file, sep='\t', index=False,
                                                                       header=False)
        command = f"bedtools intersect -wb -a {self.annotations} -b {self.temp_bed_file} " \
                  f"|cut -f1,3,5,9,13 > {self.temp_bed_results_file}"
        self.logger.info(f'bedtools command: {command}')
        return_code = self.check_call(command)
        if return_code != 0:
            self.logger.error(f"Failed to run bedtools. "
                            f"Return code: {return_code}")
            sys.exit(1)
        else:
            self.logger.info('Finished finding overlaps between barcodes and features')
            # todo leave this file, it turned out to be useful for SATAY
            #  -> if want to summarize over different intervals
            # os.remove(self.temp_bed_file)

    def _find_closest_feature(self):
        """
        """
        assert self.annotations # todo does this actually do anything?
        self.logger.info('------------------')
        self.logger.info('Annotating mapped barcodes. Finding closest feature for each barcode')
        bed_map = self.positions.copy().reset_index()
        bed_map['startOffBy1'] = bed_map['sstart'] - 1
        bed_map[['sseqid', 'startOffBy1', 'sstart', 'barcode']].to_csv(self.temp_bed_file, sep='\t', index=False,
                                                                       header=False)
        tmp_annotations = self.annotations.with_suffix('.sorted.gff')
        tmp_bed = self.temp_bed_file.with_suffix('.sorted.bed')
        command = f"bedtools sort -i {self.temp_bed_file} > {tmp_bed};" \
                  f"bedtools sort -i {self.annotations} |grep ID=gene > {tmp_annotations};" \
                  f"bedtools closest -a {tmp_annotations} -b {tmp_bed} -d " \
                  f"|cut -f1,3,5,9,13,14 > {self.temp_bed_closest_file}"
        self.logger.info(f'Bedtools command: {command}')
        return_code = self.check_call(command)
        if return_code != 0:
            self.logger.error(f"Failed to run bedtools. "
                              f"Return code: {return_code}")
            sys.exit(1)
        else:
            self.logger.info('Finished finding features close to or overlapping barcodes')
            os.remove(tmp_bed)
            os.remove(tmp_annotations)
            # os.remove(self.temp_bed_file)

    def _add_bedintersect_results_to_positions(self, feature_type, intersect=True):
        """
        takes output_map file produced by _find_annotation_overlaps and merges it with barcode map on barcode
        """
        assert self.identifiers is not None
        if intersect:
            antd_positions = pd.read_table(self.temp_bed_results_file, header=None)
            antd_positions.columns = ['seq', 'feature', 'position', 'gene_info', 'barcode']
            antd_positions = antd_positions[antd_positions.feature == feature_type]
            antd_positions['distance_to_feature'] = 0
        else:
            antd_positions = pd.read_table(self.temp_bed_closest_file, header=None)
            antd_positions.columns = ['seq', 'feature', 'position', 'gene_info', 'barcode', 'distance_to_feature']
            antd_positions = antd_positions[antd_positions.feature == feature_type]

            def get_all_genes_with_min_distance(df):
                minv = df['distance_to_feature'].min()
                return df[df['distance_to_feature'] == minv]

            antd_positions = (antd_positions.groupby('barcode')
                              .apply(get_all_genes_with_min_distance)
                              .drop(['barcode'], axis=1).reset_index())

        for feat in self.identifiers:
            pattern = f'({feat}=.+?;|{feat}=.+?$)'
            antd_positions[feat] = (antd_positions['gene_info']
                                    .str.extract(pattern, expand=False)
                                    .str.replace(f'{feat}=', '')
                                    .str.strip(';'))
            if antd_positions[feat].isna().all():
                self.logger.warning(f'"{feat}" identifier not found')
        antd_positions = antd_positions[list(self.identifiers) + ['barcode', 'distance_to_feature']]
        self.annotated_positions = (self.positions.copy()
                                    .merge(antd_positions, how='left', on='barcode'))
        # os.remove(self.temp_bed_results_file)

    def map_insertions(self, min_host_bases=20, filter_below=100):
        self.extract_barcodes(min_host_bases=min_host_bases)
        self._write_barcodes_to_fasta()
        self._blast_barcode_host(4)
        self._find_most_likely_positions(filter_below=filter_below)
        self._merge_colliding_barcodes()

    def annotate(self, annotations_file, intersect=True, feature_type='gene', identifiers=('Name', 'locus_tag')):
        assert not self.positions.empty
        self.annotations = Path(annotations_file)
        self.identifiers = identifiers
        if intersect:
            self._find_annotation_overlaps()
        else:
            self._find_closest_feature()
        self._add_bedintersect_results_to_positions(feature_type, intersect)

    def write_mapfile(self):
        rename_dict = {'sstart': 'insertion_site',
                       'sseqid': 'chr',
                       'sstrand': 'strand',
                       'total_count': 'number_of_reads'}
        self.positions = self.positions.rename(rename_dict, axis=1)
        self.positions.to_csv(self.map_file)
        if not self.annotated_positions.empty:
            self.annotated_positions = self.annotated_positions.rename(rename_dict, axis=1)
            self.annotated_positions.to_csv(self.map_file.with_suffix('.annotated.csv'), index=False)
