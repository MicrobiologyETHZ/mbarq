from tnseq2.src.sequence import stream_fa, FastA
from tnseq2.src.extract_barcodes import extract_barcode_host
from tnseq2.src.commandline import check_call
from typing import List, Tuple, Generator,  Iterator
import logging
import collections
import sys
import pandas as pd
from pathlib import Path
from logging import Logger
from tnseq2.src.tnseq_logger import get_logger


def shutdown(status=0):
    '''
    Generic shutdown commands
    '''
    logging.info('Finishing tnseq pipeline with status:\t{}'.format(status))
    sys.exit(status)


def extract_barcodes(inserts: Generator[Tuple[FastA, FastA], None, None],
                    tp2: str='GTGTATAAGAGACAG', bc2tp2:int=13, bcLen:int=17, before:bool=True,
                      min_host_bases: int = 20, logger: Logger = logging.getLogger()) -> List[Tuple[str, str]]:
    """
    -----|BARCODE|----------|TN end sequence (tp2)|---Host------
    -----|-bcLen-|--bc2tp2--|---------tp2---------|-------------
    -----|-17bp--|---13bp---|---------15bp--------|----?--------
    ---(-30)---(-13)-------(0)---------------------------------
    :param inserts genreator over r1 and r2 sequence.FastA

    """

    if not tp2:
        logging.error('Unkown transposon')
        shutdown(1)

    total_inserts: int = 0
    inserts_with_tp2: int = 0
    inserts_without_tp2: int = 0
    inserts_with_tp2_but_short_bc: int = 0
    inserts_with_tp2_with_good_bc: int = 0
    inserts_with_tp2_with_good_bc_short_host: int = 0
    sequences: List[Tuple[str, str]] = []
    r1: FastA = None
    r2: FastA = None

    for total_inserts, (r1, r2) in enumerate(inserts):
        if total_inserts % 100000 == 0:
            logging.info(f'Processed {total_inserts} reads')
        if tp2 in r1.sequence:
            inserts_with_tp2 += 1
            barcode, host_sequence = extract_barcode_host(r1, tp2, bc2tp2, bcLen, before)
            if barcode:
                inserts_with_tp2_with_good_bc += 1
                if len(host_sequence) < min_host_bases:
                    inserts_with_tp2_with_good_bc_short_host += 1
                else:
                    sequences.append((barcode, host_sequence))
            else:
                inserts_with_tp2_but_short_bc += 1
        else:
            inserts_without_tp2 += 1
    logger.info(f'Processed {total_inserts} reads')
    logger.info('Extraction statistics:')
    logger.info(f'\tTotal inserts:\t{total_inserts}\t100.0%')
    logger.info(f'\tInserts w transposon:\t{inserts_with_tp2}\t{int(100 * ((inserts_with_tp2 * 100.0)/total_inserts))/100.0}%')
    logger.info(f'\t\tand w barcode:\t{inserts_with_tp2_with_good_bc}\t{int(100 * ((inserts_with_tp2_with_good_bc * 100.0)/total_inserts))/100.0}%')
    logger.info(f'\t\t\tand good host sequence:\t{len(sequences)}\t{int(100 * ((len(sequences) * 100.0)/total_inserts))/100.0}% --> Used for downstream analysis')
    logger.info(f'\t\t\tand short host sequence:\t{inserts_with_tp2_with_good_bc_short_host}\t{int(100 * ((inserts_with_tp2_with_good_bc_short_host * 100.0)/total_inserts))/100.0}%')
    logger.info(f'\t\tand w/o barcode:\t{inserts_with_tp2_but_short_bc}\t{int(100 * ((inserts_with_tp2_but_short_bc * 100.0)/total_inserts))/100.0}%')
    logger.info(f'\tInserts w/o transposon:\t{inserts_without_tp2}\t{int(100 * ((inserts_without_tp2 * 100.0)/total_inserts))/100.0}%')
    return sequences


def prepare_extract_barcodes(in_r1_file: str, in_r2_file: str, temp_fasta_file: str,
                             tp2: str='GTGTATAAGAGACAG', bc2tp2:int=13, bcLen:int=17, before:bool=True,
                             min_host_bases: int = 20, logger: Logger = logging.getLogger()) -> None:

    """
    Wrapper function to extract barcodes and host sequences from
    sequence file. Sequences will then be written to fasta file
    for downstream analysis
    :param in_r1_file:
    :param in_r2_file:
    :param temp_fasta_file:
    :return:
    """

    fq1_stream: Generator[FastA, None, None] = stream_fa(in_r1_file)
    fq2_stream: Generator[FastA, None, None] = stream_fa(in_r2_file)
    logger.info('----------------')
    logger.info('Step 1.1: Start extraction')
    inserts: Iterator[Tuple[FastA, FastA]] = zip(fq1_stream, fq2_stream)
    alignments = extract_barcodes(inserts, tp2, bc2tp2, bcLen, before, min_host_bases, logger)
    logger.info('Step 1.1: Finished extraction')
    logger.info('----------------')
    logger.info('Step 1.2: Start writing dereplicated barcode/host pairs.')
    with open(temp_fasta_file, 'w') as handle:
        barcode_2_sequences = collections.defaultdict(list)
        for alignment in alignments:
            barcode_2_sequences[alignment[0]].append(alignment[1])
        tot = 1
        for barcode, sequences in barcode_2_sequences.items():
            for sequence, cnt in collections.Counter(sequences).most_common():
                handle.write(f'>{tot}_bc_{barcode}_cnt_{cnt}\n{sequence}\n')
                tot += 1
    logger.info(f'Step 1.2: Finished writing {tot-1} dereplicated barcode/host pairs.')
    logger.info('----------------')


def map_host_map(temp_fasta_file: str, temp_blastn_file: str, genome: str = '',
                 blastdb: str = '', blast_threads: int = 1, logger: Logger = logging.getLogger()) -> Tuple[int, int]:
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
    if not blastdb and not genome:
        logging.error("Genome or blast database needs to be provided")
        sys.exit(1)
    elif not blastdb:
        command0 = f'makeblastdb -in {genome} -dbtype nucl'
        logging.info('No blastdb provided, generating blast db from genome.')
        logging.info(f"Blast command:\t{command0}")
        db_return_code = check_call(command0)
        blastdb = genome
    command = f'blastn -task blastn -db {blastdb} -out {temp_blastn_file} -query {temp_fasta_file} ' \
              f'-outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sstrand" ' \
              f'-num_threads {blast_threads}'
    logger.info(f'Blastn command:\t{command}')
    blast_return_code = check_call(command)
    return db_return_code, blast_return_code


def new_map_annotate(blast_file: str, filter_below: int = 100,  logger: Logger = logging.getLogger()) -> pd.DataFrame:
    """
    Takes in blast file, and provides most likely locations for each barcode
    :param: blast_file
    :param: filter_below
    :param: logger
    :return: pd.DataFrame
    """
    df = pd.read_table(blast_file, header=None)
    df.columns = "qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sstrand".split()
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
    best_hits = best_hits[best_hits.cnt > filter_below]  # todo should be 100 or optional
    logger.info(f"Number of barcodes: {best_hits.barcode.nunique()}")
    # Create best hits data frame by merging best_hits with other columns from blast file
    # There still could be multiple hits for each qseqid, if they have the same blast score
    queryBH = best_hits.merge(df, how='left', on=['qseqid', 'bitscore'])
    # Check for multimapping
    multimap = (queryBH.groupby(['barcode']).sstart.std(ddof=0) > 5).reset_index().rename({'sstart': 'multimap'}, axis=1)
    queryBH = queryBH.merge(multimap, on='barcode')
    # For each barcode select the position supported by most reads
    queryBH = queryBH.sort_values(['barcode', 'cnt'], ascending=False)
    queryBH['rank'] = queryBH.groupby(['barcode']).cumcount()
    queryBH = queryBH[queryBH['rank'] == 0].copy()
    queryBH.drop('rank', axis=1, inplace=True)
    return queryBH


def editdistance(seq1: str, seq2: str) -> int:
    '''
    Calculate the edit distance between 2 sequences with identical length.
    Will throw an error if the length of both sequences differs
    :param seq1:
    :param seq2:
    :return:
    '''

    if seq1 == seq2:
        return 0
    if len(seq1) != len(seq2):
        raise Exception(f'{seq1} and {seq2} have different length. Edit distance can be computed on same length sequences only.')
    dist = 0
    for letter1, letter2 in zip(seq1, seq2):
        if letter1 != letter2:
            dist += 1
    return dist


def merge_colliding_bcs(potential_positions: pd.DataFrame, edit_dist=3, logger: Logger = logging.getLogger()) -> pd.DataFrame:
    """

    Takes output of new_map_annotate, and merges colliding barcodes

    """
    pps = potential_positions[['sseqid', 'sstart', 'send', 'sstrand', 'barcode', 'cnt', 'total_count', 'multimap']]
    # Find barcodes mapped to the same position
    collisions = (potential_positions.groupby(['sseqid', 'sstart']).barcode.nunique()
                  .reset_index().rename({'barcode': 'num_uniq_bc'}, axis=1))
    # Get all the unique barcodes mapped to unique positions
    final_bcs = list(collisions[collisions.num_uniq_bc == 1].merge(pps, on=['sseqid', 'sstart']).barcode.values)
    # Collisions are only positions with more than one bc mapped
    collisions = collisions[collisions.num_uniq_bc > 1]
    collisions = collisions.merge(pps, on=['sseqid', 'sstart'])
    logger.info(f'Number of collissions: {collisions.shape[0]}')
    for position, df in collisions.groupby(['sseqid', 'sstart']):  # todo refactor
        # Get the barcode that was seen the most
        bc_idx = df.total_count.astype(int).idxmax()
        # That is likely the correct barcode
        likely_barcode = df.loc[bc_idx].barcode
        all_barcodes = list(df.barcode.values)
        all_barcodes.remove(likely_barcode)
        final_bcs.append(likely_barcode)

        # Check if edit distance of all other barcodes is > than specified, if yes, add it to the list as well
        # (or should ignore it? Different list? is this where other library contamination comes in?)
        for bc in all_barcodes:
            ced = editdistance(bc, likely_barcode)
            if ced > edit_dist:
                final_bcs.append(bc)
    barcode_map = pps[pps.barcode.isin(final_bcs)][
        ['barcode', 'cnt', 'sstart', 'send', 'sseqid', 'sstrand', 'multimap']]
    barcode_map = barcode_map.set_index('barcode')
    return barcode_map


def add_gff_annotations(barcode_map: pd.DataFrame, bed_file: str, gff_file: str, output_map: str) -> int:

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
    bed_map = barcode_map.copy().reset_index()
    bed_map['startOffBy1'] = bed_map['sstart'] - 1
    # todo make creation of bedfile internal, remove after done
    bed_map[['sseqid', 'startOffBy1', 'startOffBy1', 'barcode']].to_csv(bed_file, sep='\t', index=False, header=False)
    command = f"bedtools intersect -wb -a {gff_file} -b {bed_file} |cut -f1,5,9,13|grep 'ID=gene' > {output_map}"
    bed_return_code = check_call(command)
    return bed_return_code


def add_gene_annotations(bedfile_path, barcode_map, barcode_file):
    '''
    takes output_map file produced by add_gff_annotations and merges it with barcode map on barcode
    '''
    bedfile = pd.read_table(bedfile_path, header=None)
    bedfile.columns = ['seq', 'position', 'gene_info', 'barcode']
    bedfile['ShortName'] = bedfile['gene_info'].str.extract(r'(Name=.+?;)', expand=False).str.replace('Name=', '').str.strip(';')
    bedfile['locus_tag'] = bedfile['gene_info'].str.extract(r'(locus_tag=.+?$)', expand=False).str.replace(
        'locus_tag=', '').str.strip(';')
    potential_positions = barcode_map.copy()
    final_map = Path(barcode_file).with_suffix('.annotated.csv')
    fdf = potential_positions.merge(bedfile[['barcode', 'ShortName', 'locus_tag']], how='left', on='barcode')
    fdf.to_csv(final_map, index=False)


def map_host_new(temp_fasta_file: str, temp_blastn_file: str, genome: str,
                 blast_threads: int, temp_bed_file: str, output_map: str, barcode_file: str, gff_file: str='', filter_below: int=100,
                 logger: Logger = logging.getLogger()) -> None:
    """
    Takes fasta file produced by prepare_extract_barcodes, runs blast and annotates

    :param temp_fasta_file: fasta file producted by prepare_extract_barcodes
    :param temp_blastn_file: file to write blast results to
    :param genome: path to the refernce genome
    :param blast_threads:
    :param temp_bed_file:
    :param output_map:
    :param barcode_file:
    :param gff_file:
    :param filter_below:
    :param logger:
    :return:
    """

    logger.info('----------------')
    logger.info('Step 2.1: Start blastn alignment of dereplicated barcode/host pairs')
    # Runs blast
    map_host_map(temp_fasta_file, temp_blastn_file, genome=genome, blast_threads=blast_threads, logger=logger)
    logger.info('Step 2.1: Finished blastn alignment of dereplicated barcode/host pairs')
    logger.info('----------------')
    logger.info('Step 2.2: Start annotating barcodes')
    # Figure out the positions
    potential_positions = new_map_annotate(temp_blastn_file, filter_below, logger)
    logger.info('Resolving barcode collisions')
    final_positions = merge_colliding_bcs(potential_positions, logger=logger)
    final_positions.to_csv(barcode_file)
    logger.info(f'GFF file provided: {gff_file}')
    if gff_file:
        logger.info('Adding gene annotations from the gff file')
        add_gff_annotations(final_positions, temp_bed_file,  gff_file,  output_map)
        add_gene_annotations(output_map, final_positions, barcode_file)
    logger.info('Finished Annotating Barcodes')
    return None


def map(r1_file, r2_file, name, output_dir, transposon, genome, gff_file,  min_host_bases=20, blast_threads=1, filter_below=100):
    logger = get_logger('Mapping_Logger', Path(output_dir) / 'tnseq2_mapping.log')
    logger.info(f'Transposon: {transposon}')
    logger.info(f"FASTQ file: {r1_file}")
    logger.info(f"GFF file: {gff_file}")
    logger.info(f"Reference Genome: {genome}")
    fasta_file = Path(output_dir)/f'{name}.fasta'
    blastn_file = Path(output_dir)/f'{name}.blastn'
    tp2 = transposon.split(':')[0]
    bc2tp2 = int(transposon.split(':')[2])
    bcLen = int(transposon.split(':')[1])
    before = True if transposon.split(':')[3] == 'before' else False
    barcode_file = Path(output_dir)/f'{name}.barcode_map.csv'
    bed_file = Path(output_dir)/f'{name}.temp.bed'
    output_map = Path(output_dir)/f'{name}.output.bed'
    logger.info("Starting Barcode Mapping")
    logger.info('Extracting Barcodes')
    prepare_extract_barcodes(r1_file,  r2_file, fasta_file,  tp2, bc2tp2, bcLen, before, min_host_bases, logger)
    map_host_new(fasta_file, blastn_file, genome, blast_threads, bed_file, output_map,  barcode_file, gff_file, filter_below, logger)





if __name__ == '__main__':

    r1 = sys.argv[1]
    r2 = sys.argv[2]
    name = sys.argv[3]
    output_dir = sys.argv[4]
    genome = sys.argv[5]
    gff_file = sys.argv[6]
    blast_threads = int(sys.argv[7])
    transposon = sys.argv[8]
    #min_host_bases = int(sys.argv[])

    map(r1, r2, name, output_dir, transposon, genome, gff_file, blast_threads)

