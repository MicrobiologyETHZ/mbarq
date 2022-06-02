import subprocess
import shlex
from pathlib import Path
import click
import sys
import logging
from mbarq.mapper import Mapper, AnnotatedMap
from mbarq.counter import BarcodeCounter
from mbarq.analysis import CountDataSet, Experiment
import sys

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
@click.option('--gff', '-a', default='', help='annotation file in GFF format', metavar='FILE')
@click.option('--name', '-n', default='', help='unique library name, '
                                               'by default will try to use FASTQ filename', metavar='STR')
@click.option('--transposon', '-tn',
              default="B17N13GTGTATAAGAGACAG",
              help="""\b
                     transposon construct structure, consisting of the following:
                       1. barcode length, written as B[# of nt], eg. B17
                       2. conserved sequence motif, usually part of transposons inverted repeat (IR), eg. GTGTATAAGAGACAG
                       3. if there are extra nucleotides between barcode and 
                          conserved sequence motif, indicate with N[# of nt], eg. N13
                    The default represents the following construct:
                    ----------------------------------------------------------------------------
                    Read      ||AGTACTTTACTACTACT||TACCTGACCGTAA||GTGTATAAGAGACAG||TTACCTGACCGAC
                    ----------||-----------------||-------------||---------------||-------------
                    Components||     barcode     ||   spacer    ||    conserved  ||     host    
                              ||                 ||             ||    motif (IR) ||               
                    ----------||-----------------||-------------||---------------||-------------
                    Encoding  ||       B17       ||     N13     ||GTGTATAAGAGACAG||              
                    ----------------------------------------------------------------------------
                    Note: relative position of barcode and conserved sequence motif matters, i.e. if conserved sequence motif comes before the barcode,
                    it should be written as GTGTATAAGAGACAGN13B17.
                     [B17N13GTGTATAAGAGACAG]
              
                   """, metavar='STR')
@click.option('--out_dir', '-o', default='.', help='output directory [.]', metavar="DIR")
@click.option('--filter_low_counts', '-l', default=0,
              help='filter out barcodes supported by [INT] or less reads [0]', metavar="INT")
#@click.option('--blast_threads', '-t', default=4, help='Blast Threads')
@click.option('--feat_type', '-ft', default='gene',
              help='feature type in the GFF file to be used for annotation, e.g. gene, exon, CDS [gene]',
              metavar="STR")
@click.option('--attributes', default='ID,Name,locus_tag',
              help='Feature attributes to extract from GFF file [ID,Name,locus_tag]', metavar="STR[,STR]")
@click.option('--closest_gene',  is_flag=True,
              help='for barcodes not directly overlapping a feature, report the closest feature [False]')
def map(forward, gff, name, transposon, out_dir, genome, filter_low_counts,
        feat_type, attributes, closest_gene):
    identifiers = tuple(attributes.split(','))
    mapper = Mapper(forward, transposon, genome=genome, name=name, output_dir=out_dir)
    mapper.map_insertions(filter_below=filter_low_counts)
    if gff:
        annotated_map = AnnotatedMap(map_file=mapper.map_file, annotation_file=gff,
                                     feature_type=feat_type,
                                     identifiers=identifiers, output_dir=out_dir,
                                     name=name)
        annotated_map.annotate(intersect=not closest_gene)


#####################
#  ANNOTATE MAPPED  #
#####################

@main.command(cls=DefaultHelp, short_help="\tannotate mapped barcodes", options_metavar='<options>')
@click.option('--barcode_file', '-i',  help='unannotated barcode file produced by "map"', metavar='FILE')
@click.option('--gff', '-a',  help='annotation file in gff format', metavar='FILE')
@click.option('--name', '-n', default='', help='unique library name, '
                                               'by default will try to use FASTQ filename', metavar='STR')
@click.option('--out_dir', '-o', default='.', help='output directory [.]', metavar="DIR")
@click.option('--feat_type', '-ft', default='gene',
              help='feature type in the gff file to be used for annotation, e.g. gene, exon, CDS [gene]',
              metavar="STR")
@click.option('--attributes',  default='ID,Name,locus_tag',
              help='Feature attributes to extract from annotation file [ID,Name,locus_tag]', metavar="STR[,STR]")
@click.option('--closest_gene', is_flag=True,
              help='for barcodes not directly overlapping a feature, report the closest feature [False]')
def annotate_mapped(barcode_file, gff, name,  out_dir, feat_type, attributes, closest_gene):
    identifiers = tuple(attributes.split(','))
    annotatedMap = AnnotatedMap(map_file=barcode_file, annotation_file=gff,
                                feature_type=feat_type, identifiers=identifiers,
                                name=name, output_dir=out_dir)
    annotatedMap.annotate(intersect=not closest_gene)

###########
#  COUNT  #
###########
@main.command(cls=DefaultHelp, short_help="\tcount barcodes in a sequencing file", options_metavar='<options>')
@click.option('--forward', '-f', required=True,
              help='input file for reads in forward orientation; FASTQ formatted; gz is ok.',
              metavar='FILE')
@click.option('--mapping_file', '-m', default='',
              help='Mapping file produced by `mbarq map`. Alternatively, will accept any csv file, where  '
                   'the first column is titled "barcode", and contains the barcode sequences. \n\n'
                   'Example: \n\n'
                   'barcode,barcodeID\n'
                   'AGACCAGTACATGACGGGTATCTCTCTGCCACTCCTGTAT,Tag_1\n '
                   '[optional] ', metavar='FILE')
@click.option('--out_dir', '-o', default='.', help='output directory [.]', metavar="DIR")
@click.option('--name', '-n', default='', help='unique library name, '
                                               'by default will try to use FASTQ filename', metavar='STR')
@click.option('--transposon', '-tn',
              default="B17N13GTGTATAAGAGACAG",
              help="""\b
                     transposon construct structure, consisting of the following:
                       1. barcode length, written as B[# of nt], eg. B17
                       2. conserved sequence motif, usually part of transposons inverted repeat (IR), eg. GTGTATAAGAGACAG
                       3. if there are extra nucleotides between barcode and 
                          conserved sequence motif, indicate with N[# of nt], eg. N13
                    The default represents the following construct:
                    ----------------------------------------------------------------------------
                    Read      ||AGTACTTTACTACTACT||TACCTGACCGTAA||GTGTATAAGAGACAG||TTACCTGACCGAC
                    ----------||-----------------||-------------||---------------||-------------
                    Components||     barcode     ||   spacer    ||    conserved  ||     host    
                              ||                 ||             ||    motif (IR) ||               
                    ----------||-----------------||-------------||---------------||-------------
                    Encoding  ||       B17       ||     N13     ||GTGTATAAGAGACAG||              
                    ----------------------------------------------------------------------------
                    Note: relative position of barcode and conserved sequence motif matters, i.e. if conserved sequence motif comes before the barcode,
                    it should be written as GTGTATAAGAGACAGN13B17. For WISH data, use GGAGGTTCACAATGTGGGAGGTCAB40
                     [B17N13GTGTATAAGAGACAG]
                     
                    """, metavar='STR')
@click.option('--edit_distance', '-e', default=2,
              help='merge barcodes with edit distances <= [INT] [2]', metavar="INT")
def count(forward, mapping_file, out_dir, transposon, name, edit_distance):
    counter = BarcodeCounter(forward, transposon, name=name, mapping_file=mapping_file,
                             output_dir=out_dir, edit_distance=edit_distance)
    counter.count_barcodes()


###########
#  MERGE  #
###########
@main.command(cls=DefaultHelp, short_help="\tmerge count files", options_metavar='<options>')
@click.option('--input_files', '-i', default='', help='list of mBARq count files to merge (comma separated)', metavar='FILE[,FILE]')
@click.option('--count_dir', '-d', default='',  help='merge all files in the directory', metavar='DIR')
@click.option('--out_dir', '-o', default='.', help='output directory', metavar='DIR')
@click.option('--name', '-n',  help='output file prefix', metavar='STR')
@click.option('--attribute', '-a',  default='',
              help='Feature attribute to keep in the merged file (ex. ID, Name, locus_tag)', metavar="STR")
def merge(input_files, count_dir, name, attribute, out_dir):
    if not input_files and not count_dir:
        sys.exit('Please specify files to merge or a directory with count files')
    elif input_files:
        input_files = [Path(f) for f in input_files.split(',')]
    else:
        input_files = [f for f in Path(count_dir).iterdir() if 'mbarq_counts' in f.name]
    count_dataset = CountDataSet(count_files=input_files, name=name,
                                 gene_column_name=attribute, output_dir=out_dir)
    count_dataset.create_count_table()


@main.command(cls=DefaultHelp, short_help="analyze transposons for differential abundance. Under construction",  options_metavar='<options>')
@click.option('--count_file', '-i',  help='count file produced by "merge"', metavar='FILE')
@click.option('--sample_data', '-s',  help='sample data', metavar='FILE')
@click.option('--control_file', '-c',  help='control barcode file', metavar='FILE')
@click.option('--gene_name', '-g',  default='Name', help='Name of the column containing gene '
                                                         'identifiers in the count file', metavar='STR')
@click.option('--treatment_column',  help='column in sample data indicating treatmen', metavar='STR')
@click.option('--batch_column',  help='column in sample data indicating batch', metavar='STR')
@click.option('--baseline',  help='treatment to use as control/baseline, ex. day 0', metavar='STR')
@click.option('--name', '-n', default='', help='experiment name, '
                                               'by default will try to use count file name', metavar='STR')
@click.option('--out_dir', '-o', default='.', help='output directory', metavar='DIR')
def analyze(count_file, sample_data, gene_name, control_file, name,
            treatment_column, baseline, batch_column, out_dir):
    exp = Experiment(count_file, sample_data, control_file, name, gene_name, treatment_column,
                     baseline, batch_column, 0.8, out_dir)
    exp.run_experiment()


if __name__ == "__main__":
    main()
