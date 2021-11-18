import subprocess
import shlex
from pathlib import Path
import click
import sys
import logging
from tnseq2.src.demultipex import demux_tnseq
from tnseq2.src.quantify import quantify
from tnseq2.src.mapping import map
from tnseq2.src.merge_counts import final_merge
@click.group()
def main():
    pass


#DEMUX
@main.command(help='demultiplex RBSeq fastq files')
@click.option('--config', '-c', help='Configuration File')
@click.option('--input_file', '-i',  help='Input FASTQ to demultiplex')
@click.option('--demux_file', '-d', help='Barcode Map, tab delimited, ex.\n\nACCT\tSample1\n\nAAGG\tSample2\n')
@click.option('--out_dir', '-o', default='.', help='Output Directory')
@click.option('--transposon', '-tn', default="GTGTATAAGAGACAG:17:13:before", help='Construct Structure:\n\n'
                                                                                'TN sequence:BC length:length of spacer between BC and TN sequence:BC position relative to TN \n\n'
                                                                                  'Default: GTGTATAAGAGACAG:17:13:before')
@click.option('--name', '-n',  default='', help="Sample Name")
@click.option('--rc', is_flag=True, help="Reverse complement the barcodes")
@click.option('--dry',  is_flag=True, help="Show commands without running them")
@click.option('--local',  is_flag=True, help="Run on local machine")
def demux(config, input_file, demux_file, out_dir, rc, transposon, dry, local, name):
    if (config or input_file) and not (config and input_file):
        if config:
            print(f'Your provided a config file: {config}')
            #click.echo("Running {}".format('locally' if local else ('dry' if dry else 'on cluster')))
            #cmd = snakemake_cmd(config, 'demux_all', dry, local)
            #click.echo(" ".join(cmd))
        else:
            print(f"You've provided a FASTQ file: {input_file}")
            demux_tnseq(input_file, demux_file, out_dir, name, transposon, rc)
    else:
        print('Provide either config or FASTQ file, not both')
        sys.exit(1)


# MAPPING
@main.command(short_help="map barcodes to reference genome")
@click.option('--forward', '-f',  help='Forward Reads')
@click.option('--reverse', '-r',  help='Reverse Reads')
@click.option('--gff', '-a', default='', help='Annotation File in gff format. '
                                              'The tool will extract Name and locus_tag from gene features')
@click.option('--genome', '-g',  help='Reference genome (FASTA)')
@click.option('--name', '-n', help='Unique library name')
@click.option('--out_dir', '-o', default='.', help='Output directory')
@click.option('--filter_low_counts', '-l', default=100, help='Filter out barcodes supported by [int] or less reads. Default: [100]')
@click.option('--blast_threads', '-t', default=1, help='Blast Threads')
@click.option('--transposon', '-tn', default="GTGTATAAGAGACAG:17:13:before", help='Construct Structure:\n\n'
                                                                                'TN sequence:BC length:length of spacer between BC and TN sequence:BC position relative to TN \n\n'
                                                                                  'Default: GTGTATAAGAGACAG:17:13:before')
def maplib(forward, reverse, gff, name, transposon, out_dir, genome, blast_threads, filter_low_counts):
    if not name:
        name = Path(forward.strip('.gz')).stem
    map(forward, reverse, name, out_dir, transposon, genome, gff_file=gff,
        blast_threads=blast_threads, filter_below=filter_low_counts)


# COUNT
@main.command(short_help="count barcodes in one sample")
@click.option('--forward', '-f',  help='Input FASTQ to count. (Forward Reads only)')
@click.option('--mapping_file', '-m', help='Barcode Map in csv format. First column must be titled "barcode", and contain the barcodes. \n\n'
                                           'Example: \n\n'
                                           'barcode,barcodeID\n\n'
                                           'AGACCAGTACATGACGGGTATCTCTCTGCCACTCCTGTAT,Tag_1\n\n')

@click.option('--out_dir', '-o', default='.', help='Output Directory')
@click.option('--sample_name', '-n', default='', help='Sample Name')
@click.option('--transposon', '-tn', default="GTGTATAAGAGACAG:17:13:before", help='Construct Structure:\n\n'
                                                                                'TN sequence:BC length:length of spacer between BC and TN sequence:BC position relative to TN (before or after) \n\n'
                                                                                  '-|BARCODE|-spacer-|--TN sequence--|-\n\n'
                                                                                  '-|-17bp--|--13bp--|GTGTATAAGAGACAG|-\n\n'
                                                                                  'Default: GTGTATAAGAGACAG:17:13:before')
def count(forward, mapping_file, out_dir, transposon, sample_name):
    if not sample_name:
        sample_name = Path(forward.strip('.gz')).stem
    quantify(forward, transposon, mapping_file, out_dir, sample_name)



# # Custom
# @main.command()
# @click.option('--config', '-c', default='configs/map_config.yaml', help='Configuration File')
# @click.option('--local',  is_flag=True, help="Run on local machine")
# @click.option('--dry',  is_flag=True, help="Show commands without running them")
# @click.option('--method', '-m',  help='Run custom command, for testing mode')
# def custom(config, local, dry, method):
#     click.echo("Mapping Barcode Libraries")
#     click.echo(f"Config file: {config}")
#     click.echo("Samples found: ")
#     click.echo("Running {}".format('locally' if local else 'on cluster'))
#     cmd = snakemake_cmd(config, method, dry, local)
#     click.echo(" ".join(cmd))
#
#


#MERGE
@main.command( short_help ="merge counts from multiple samples. Under construction.")
@click.option('--count_dir', '-d',  help='Input directory with count files')
@click.option('--meta_file', '-m', default='', help='Meta file, format: ...')
@click.option('--control_file', '-b', default='', help="File with WITS info, format: ...")
#@click.option('--out_file', '-o', default='', help='Output Directory')
@click.option('--runid', '-n', default='', help='Run/Experiment Name')
def merge(config, count_dir, meta_file, control_file, runid):
    logging.info(f"You've provided a counts directory: {count_dir}")
    final_merge(count_dir,  meta_file, control_file, runid)



@main.command(short_help="analyze transposons for differential abundance. Under construction")
@click.option('--config', '-c', default='configs/analyze_config.yaml', help='Configuration File')
def analyze():
    pass

if __name__ == "__main__":
    main()
