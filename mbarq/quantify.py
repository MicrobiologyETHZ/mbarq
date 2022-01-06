from tnseq2.src.sequence import stream_fa, FastA
from tnseq2.src.extract_barcodes import extract_barcode_host
from tnseq2.src.tnseq_logger import get_logger
from tnseq2.src.commandline import *
from typing import List, Tuple, Generator, Dict, Iterator, DefaultDict, Counter, Set
import logging
import collections
import sys
import pandas as pd
from pathlib import Path
from logging import Logger


def quantify_load_fq_barcodes(in_fq_file: str, tp2: str = 'GTGTATAAGAGACAG',
                              bc2tp2: int = 13, bcLen: int = 17, before: bool = True,
                              logger: Logger = logging.getLogger()) -> Counter[str]:
    '''
    Count barcodes in a sequencing file (fasta/fastq)
    This currently the default:

    -----|BARCODE|----------|TN end sequence (tp2)|---Host------
    -----|-bcLen-|--bc2tp2--|---------tp2---------|-------------
    -----|-17bp--|---13bp---|---------15bp--------|----?--------
    ---(-30)---(-13)-------(0)---------------------------------

    Uses extract_barcode_host to actually get the barcode sequence

    :param in_fq_file: path to forward reads file (FASTQ/FASTA)
    :param tp2: conserved tn sequence
    :param bc2tp2: distance between barcode and conserved sequence
    :param bcLen: length of the barcode
    :param before: if the barcode before or after the conserved sequence
    :param logger:

    :return: counter[barcode]=count
    '''
    fq1_stream = stream_fa(in_fq_file)
    logger.info("Counting reads with and without transposon sequence")
    with_tp2: int = 0
    without_tp2: int = 0
    with_tp_but_short: int = 0
    cnter: Counter[str] = collections.Counter()
    for total_inserts, r1 in enumerate(fq1_stream, 1):
        if total_inserts % 1000000 == 0:
            logger.info(f'\tReads processed:\t{total_inserts}')
        if tp2 in r1.sequence:
            barcode, _ = extract_barcode_host(r1, tp2, bc2tp2, bcLen, before)
            if not barcode:
                with_tp_but_short += 1
            else:
                with_tp2 += 1
                cnter[barcode] += 1
        else:
            without_tp2 += 1
    logger.info(f'Reads processed:\t{total_inserts}')
    logger.info(f'FastA/Q Stats:')
    logger.info(f'\tTotal Reads:\t{total_inserts}')
    logger.info(f'\tReads w transposon and good barcode:\t{with_tp2}')
    logger.info(f'\tReads w transposon but short barcode:\t{with_tp_but_short}')
    logger.info(f'\tReads w/o transposon:\t{without_tp2}')
    logger.info(f'\t% Reads w/o transposon: \t {without_tp2/total_inserts*100}')
    return cnter


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
        raise Exception(
            f'{seq1} and {seq2} have different length. Edit distance can be computed on same length sequences only.')

    return sum(letter1 != letter2 for letter1, letter2 in zip(seq1, seq2))



def get_similar(bc_to_id: list, bc_mapped: list, logger: Logger):
    """

    Take list of unannotated barcodes (bc_to_id) and list mapped barcodes.
    Calculate edit distance between each unannotated barcode and each mapped barcode
    For each unannotated barcode get a mapped barcode with smallest edit distance

    How to make this faster?

    """
    # slow slow slow
    assert len(bc_to_id) > 0
    distances = {}
    logger.info(f'Calculating edit distances between annotated and unannotated barcodes')
    for total, bc in enumerate(bc_to_id):
        min_dist = 1000
        match = ''
        if total % 5000 == 0:
            logger.info(f'\tBarcodes processed:\t{total}')
        for mbc in bc_mapped:
            actual_dist = editdistance(bc, mbc)
            if actual_dist < min_dist:
                min_dist = actual_dist
                match = mbc
        distances[bc] = [min_dist, match]

    dist = pd.DataFrame(distances).T.reset_index()

    dist.columns = ['barcode', 'editdistance', 'match']
    return dist


def annotate_barcodes(cnter, barcode_map_file, logger, edit_cutoff=3):
    cnts_df = pd.DataFrame.from_dict(cnter, orient='index').reset_index()
    cnts_df.columns = ['barcode', 'barcode_cnt']
    cnts_df = cnts_df[cnts_df['barcode_cnt'] > 1]
    logger.info(f'Number of unique barcodes to annotate: {cnts_df.barcode.nunique()}')
    if cnts_df.empty:
        logger.error('No barcodes with counts > 1 found')
        sys.exit(1)
    # If no mapping file, return just the counts and an empty data frame
    if not barcode_map_file:
        logger.info('No mapping file found, skipping the annotation')
        return cnts_df, pd.DataFrame()

    bc_df = pd.read_csv(barcode_map_file)
    # This needs to be done in the mapping step. Here want to be flexible for mapping file input.
    #bc_df.columns = 'barcode,libcnt,sstart,send,sseqid,sstrand,multimap,ShortName,locus_tag'.split(',')
    to_keep = bc_df.columns
    logger.info(f'Mapping file used: {barcode_map_file}')
    logger.info(f'Columns found in the mapping file: {", ".join(bc_df.columns)}')
    filter_col = to_keep[1]
    if 'barcode' not in bc_df.columns:
        logger.error('No column "barcode" found')
        sys.exit(1)

    annotated_cnts = cnts_df.merge(bc_df, how='outer', on='barcode')
    with_ids = annotated_cnts[(annotated_cnts[filter_col].notnull()) & (annotated_cnts.barcode_cnt.notnull())]
    logger.info(f'Number of annotated barcodes: {with_ids.shape[0]}')
    no_ids = annotated_cnts[annotated_cnts[filter_col].isna()]
    logger.info(f'Number of unannotated barcodes:{no_ids.shape[0]}')
    if no_ids.shape[0] > 0:
        logger.info(f"Merging  barcodes with edit distance < {edit_cutoff}")
        distances = get_similar(no_ids.barcode.values, bc_df.barcode.values, logger)
        no_matches = distances[distances.editdistance >= edit_cutoff]
        with_matches = distances[distances.editdistance < edit_cutoff]
        with_matches = (with_matches.merge(cnts_df, how='left', on='barcode')
                        .drop(['barcode', 'editdistance'], axis=1)
                        .rename({'match': 'barcode'}, axis=1))
        never_ided = no_ids[no_ids.barcode.isin(no_matches.barcode.values)]
        logger.info(f'Number of barcodes w/o annotation: {never_ided.shape[0]}')
        all_ids = pd.concat([with_ids[['barcode', 'barcode_cnt']].drop_duplicates(), with_matches])

    else:
        logger.info(f'All barcodes annotated')
        all_ids = with_ids[['barcode', 'barcode_cnt']].drop_duplicates()
        never_ided = pd.DataFrame()
    all_ids = all_ids.groupby('barcode').barcode_cnt.sum().reset_index()
    final_ids = all_ids.merge(annotated_cnts[to_keep], on='barcode', how='left')
    if 'ShortName' in final_ids.columns:
        final_ids.ShortName.fillna(final_ids.barcode, inplace=True)
    return final_ids, never_ided.dropna(axis=1)


def quantify(fasta_file, transposon,  map_file, outdir, prefix):
    logger = get_logger('Quantify_Logger', Path(outdir)/'tnseq2_counting.log')
    logger.info(f'Transposon: {transposon}')
    tp2 = transposon.split(':')[0]
    bc2tp2 = int(transposon.split(':')[2])
    bcLen = int(transposon.split(':')[1])
    before = True if transposon.split(':')[3] == 'before' else False
    logger.info('Counting barcodes')
    cnter = quantify_load_fq_barcodes(fasta_file, tp2, bc2tp2, bcLen,  before, logger)
    logger.info('Annotating barcodes')
    cnts, no_ids = annotate_barcodes(cnter, map_file, logger, edit_cutoff=3)
    cnts.to_csv(Path(outdir) / f'{prefix}_counts_mapped.csv')
    no_ids.to_csv(Path(outdir) / f'{prefix}_counts_unmapped.csv')
    logger.info('All Done!')
    return cnts, no_ids

#
# if __name__ == "__main__":
#     #fasta_file = sys.argv[1]
#     #map_file = sys.argv[2]
#     #outdir = sys.argv[3]
#     #prefix = sys.argv[4]
#     #quantify(fasta_file, transposon, map_file, outdir, prefix)



