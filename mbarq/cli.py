import subprocess
import shlex
from pathlib import Path
import click
import sys
import logging
from mbarq.mapper import Mapper, AnnotatedMap
from mbarq.counter import BarcodeCounter


class DefaultHelp(click.Command):
    def __init__(self, *args, **kwargs):
        context_settings = kwargs.setdefault('context_settings', {})
        if 'help_option_names' not in context_settings:
            context_settings['help_option_names'] = ['-h', '--help']
        self.help_flag = context_settings['help_option_names'][0]
        super(DefaultHelp, self).__init__(*args, **kwargs)

    def parse_args(self, ctx, args):
        if not args:
            args = [self.help_flag]
        return super(DefaultHelp, self).parse_args(ctx, args)


@click.group(options_metavar='', subcommand_metavar='<command> [options]', invoke_without_command=False)
def main():
    """
    \b
    Program: mbarq - a tool for analysis of barcode mutagenesis data
    Version: 1.0.0
    Reference:

    Type mbarq <command> to print the help for a specific command

    """


#
# # DEMUX
# @main.command(help='demultiplex RBSeq fastq files')
# @click.option('--config', '-c', help='Configuration File')
# @click.option('--input_file', '-i', help='Input FASTQ to demultiplex')
# @click.option('--demux_file', '-d', help='Barcode Map, tab delimited, ex.\n\nACCT\tSample1\n\nAAGG\tSample2\n')
# @click.option('--out_dir', '-o', default='.', help='Output Directory')
# @click.option('--transposon', '-tn', default="GTGTATAAGAGACAG:17:13:before", help='Construct Structure:\n\n'
#                                                                                   'TN sequence:BC length:length of spacer between BC and TN sequence:BC position relative to TN \n\n'
#                                                                                   'Default: GTGTATAAGAGACAG:17:13:before')
# @click.option('--name', '-n', default='', help="Sample Name")
# @click.option('--rc', is_flag=True, help="Reverse complement the barcodes")
# @click.option('--dry', is_flag=True, help="Show commands without running them")
# @click.option('--local', is_flag=True, help="Run on local machine")
# def demux(config, input_file, demux_file, out_dir, rc, transposon, dry, local, name):
#     if (config or input_file) and not (config and input_file):
#         if config:
#             print(f'Your provided a config file: {config}')
#             # click.echo("Running {}".format('locally' if local else ('dry' if dry else 'on cluster')))
#             # cmd = snakemake_cmd(config, 'demux_all', dry, local)
#             # click.echo(" ".join(cmd))
#         else:
#             print(f"You've provided a FASTQ file: {input_file}")
#             demux_tnseq(input_file, demux_file, out_dir, name, transposon, rc)
#     else:
#         print('Provide either config or FASTQ file, not both')
#         sys.exit(1)


##########
#   MAP  #
##########

@main.command(cls=DefaultHelp, short_help="\tmap barcodes to reference genome", options_metavar='<options>')
@click.option('--forward', '-f', required=True,
              help='input file for reads in forward orientation; FASTQ formatted; gz is ok.',
              metavar='FILE')
@click.option('--genome', '-g', required=True, help='reference genome in FASTA format', metavar='FILE')
@click.option('--gff', '-a', default='', help='annotation file in gff format', metavar='FILE')
@click.option('--name', '-n', default='', help='unique library name, '
                                               'by default will try to use FASTQ filename', metavar='STR')
@click.option('--transposon', '-tn',
              default="B17N13GTGTATAAGAGACAG",
              help="""\b
                     transposon construct structure, consisting of the following:
                       1. conserved transposon sequence, eg. GTGTATAAGAGACAG
                       2. barcode length, written as B[# of nt], eg. B17
                       3. if there are extra nucleotides between barcode and 
                          transposon sequence, indicate with N[# of nt], eg. N13
                    Note: relative position of barcode and transpson matters, 
                    the default represents the following construct:
                    ---|BARCODE (17 nt)|--spacer (13 nt)--|GTGTATAAGAGACAG|---HOST--
                     [B17N13GTGTATAAGAGACAG]
              
                   """, metavar='STR')
@click.option('--out_dir', '-o', default='.', help='output directory [.]', metavar="DIR")
@click.option('--filter_low_counts', '-l', default=0,
              help='filter out barcodes supported by [INT] or less reads [0]', metavar="INT")
#@click.option('--blast_threads', '-t', default=4, help='Blast Threads')
@click.option('--feat_type', '-ft', default='gene',
              help='feature type in the gff file to be used for annotation, e.g. gene, exon, CDS [gene]',
              metavar="STR")
@click.option('--identifiers', '-i', default='ID,Name,locus_tag',
              help='Feature identifiers to extract from annotation file [ID,Name,locus_tag]', metavar="STR[,STR]")
@click.option('--closest_gene',  is_flag=True,
              help='for barcodes not directly overlapping a feature, report the closest feature [False]')
def map(forward, gff, name, transposon, out_dir, genome, filter_low_counts,
        feat_type, identifiers, closest_gene):
    identifiers = tuple(identifiers.split(','))
    mapper = Mapper(forward, transposon, genome=genome, name=name, output_dir=out_dir)
    mapper.map_insertions(filter_below=filter_low_counts)
    if gff:
        annotated_map = AnnotatedMap(map_file=mapper.map_file, annotation_file=gff,
                                     feature_type=feat_type,
                                     identifiers=identifiers, output_dir=out_dir)
        annotated_map.annotate(intersect=not closest_gene)


#####################
#  ANNOTATE MAPPED  #
#####################

@main.command(cls=DefaultHelp, short_help="\tannotate mapped barcodes", options_metavar='<options>')
@click.option('--barcode_file', '-i', default='', help='unannotated barcode file produced by "map"', metavar='FILE')
@click.option('--gff', '-a', default='', help='annotation file in gff format', metavar='FILE')
@click.option('--name', '-n', default='', help='unique library name, '
                                               'by default will try to use FASTQ filename', metavar='STR')
@click.option('--out_dir', '-o', default='.', help='output directory [.]', metavar="DIR")
@click.option('--feat_type', '-ft', default='gene',
              help='feature type in the gff file to be used for annotation, e.g. gene, exon, CDS [gene]',
              metavar="STR")
@click.option('--identifiers', '-i', default='ID,Name,locus_tag',
              help='Feature identifiers to extract from annotation file [ID,Name,locus_tag]', metavar="STR[,STR]")
@click.option('--closest_gene', is_flag=True,
              help='for barcodes not directly overlapping a feature, report the closest feature [False]')
def annotateMapped(barcode_file, gff, name,  out_dir, feat_type, identifiers, closest_gene):
    identifiers = tuple(identifiers.split(','))
    annotatedMap = AnnotatedMap(barcode_file, gff, feature_type=feat_type, identifiers=identifiers,
                                name=name, output_dir=out_dir)
    annotatedMap.annotate(intersect=not closest_gene)


# COUNT
@main.command(cls=DefaultHelp, short_help="\tcount barcodes in a sequencing file", options_metavar='<options>')
@click.option('--forward', '-f', required=True,
              help='input file for reads in forward orientation; FASTQ formatted; gz is ok.',
              metavar='FILE')
@click.option('--mapping_file', '-m', default='',
              help='Barcode map/annotation file in csv format.'
                   'First column must be titled "barcode", and contain the barcode sequences. [optional] \n\n'
                   'Example: \n\n'
                   'barcode,barcodeID\n'
                   'AGACCAGTACATGACGGGTATCTCTCTGCCACTCCTGTAT,Tag_1\n\n ', metavar='FILE')
@click.option('--out_dir', '-o', default='.', help='output directory [.]', metavar="DIR")
@click.option('--name', '-n', default='', help='unique library name, '
                                               'by default will try to use FASTQ filename', metavar='STR')
@click.option('--transposon', '-tn',
              default="B17N13GTGTATAAGAGACAG",
              help="""\b
                    transposon construct structure, consisting of the following:
                       1. conserved transposon sequence, eg. GTGTATAAGAGACAG
                       2. barcode length, written as B[# of nt], eg. B17
                       3. if there are extra nucleotides between barcode and 
                          transposon sequence, indicate with N[# of nt], eg. N13
                    Note: relative position of barcode and transpson matters, 
                    the default (RBSeq tn) represents the following construct:
                    ---|BARCODE (17 nt)|--spacer (13 nt)--|GTGTATAAGAGACAG|---HOST--
                     [B17N13GTGTATAAGAGACAG]
                     For WISH use GGAGGTTCACAATGTGGGAGGTCAB40
                     
                    """, metavar='STR')
@click.option('--edit_distance', '-e', default=2,
              help='merge barcodes with edit distances <= [INT] [2]', metavar="INT")
def count(forward, mapping_file, out_dir, transposon, name, edit_distance):
    counter = BarcodeCounter(forward, transposon, name=name, mapping_file=mapping_file,
                             output_dir=out_dir, edit_distance=edit_distance)
    counter.count_barcodes()


# # # Custom
# # @main.command()
# # @click.option('--config', '-c', default='configs/map_config.yaml', help='Configuration File')
# # @click.option('--local',  is_flag=True, help="Run on local machine")
# # @click.option('--dry',  is_flag=True, help="Show commands without running them")
# # @click.option('--method', '-m',  help='Run custom command, for testing mode')
# # def custom(config, local, dry, method):
# #     click.echo("Mapping Barcode Libraries")
# #     click.echo(f"Config file: {config}")
# #     click.echo("Samples found: ")
# #     click.echo("Running {}".format('locally' if local else 'on cluster'))
# #     cmd = snakemake_cmd(config, method, dry, local)
# #     click.echo(" ".join(cmd))
# #
# #
#
#
# # MERGE
# @main.command(short_help="merge counts from multiple samples. Under construction.")
# @click.option('--count_dir', '-d', help='Input directory with count files')
# @click.option('--meta_file', '-m', default='', help='Meta file, format: ...')
# @click.option('--control_file', '-b', default='', help="File with WITS info, format: ...")
# # @click.option('--out_file', '-o', default='', help='Output Directory')
# @click.option('--runid', '-n', default='', help='Run/Experiment Name')
# def merge(config, count_dir, meta_file, control_file, runid):
#     logging.info(f"You've provided a counts directory: {count_dir}")
#     final_merge(count_dir, meta_file, control_file, runid)
#
#
# @main.command(short_help="analyze transposons for differential abundance. Under construction")
# @click.option('--config', '-c', default='configs/analyze_config.yaml', help='Configuration File')
# def analyze():
#     pass


if __name__ == "__main__":
    main()
