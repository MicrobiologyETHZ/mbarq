# Analysis of large-throughput barcoded screen data

Example usage:
1. Mapping, quantification and analysis of random transposon mutagenesis (RB-Seq) experiments. 
2. Quantification of custom barcoded strains 
3. More to come! 

## Installation:

- Clone the repository and create and activate conda environment

```
git clone https://github.com/MicrobiologyETHZ/mbarq.git
cd mbarq
#conda env create -f mbarq_environment.yaml
mamba env create -f mbarq_environment.yaml
conda activate mbarq
pip install -e .
mbarq --help

```

### Testting

```shell
tar -xvzf tests/expected_outcomes.tar.gz
 tar -xvzf tests/test_files.tar.gz
 pytest tests/unit/test_mapper.py 
```

### Identify insertion sites for RB-Seq library


```shell
Usage: mbarq map <options>

Options:
  -f, --forward FILE           input file for reads in forward orientation;
                               FASTQ formatted; gz is ok.  [required]
  -g, --genome FILE            reference genome in FASTA format  [required]
  -a, --gff FILE               annotation file in gff format
  -n, --name STR               unique library name, by default will try to use
                               FASTQ filename
  -tn, --transposon STR         transposon construct structure, consisting of the following:
                                  1. conserved transposon sequence, eg. GTGTATAAGAGACAG
                                  2. barcode length, written as B[# of nt], eg. B17
                                  3. if there are extra nucleotides between barcode and 
                                     transposon sequence, indicate with N[# of nt], eg. N13
                               Note: relative position of barcode and transpson matters, 
                               the default represents the following construct:
                               ---|BARCODE (17 nt)|--spacer (13 nt)--|GTGTATAAGAGACAG|---HOST--
                                [B17N13GTGTATAAGAGACAG]
  -o, --out_dir DIR            output directory [.]
  -l, --filter_low_counts INT  filter out barcodes supported by [INT] or less
                               reads [0]
  -ft, --feat_type STR         feature type in the gff file to be used for
                               annotation, e.g. gene, exon, CDS [gene]
  -i, --identifiers STR[,STR]  Feature identifiers to extract from annotation
                               file [ID,Name,locus_tag]
  -c, --closest_gene BOOLEAN   for barcodes not directly overlapping a
                               feature, report the closest feature [False]

```
### Example usage:

Required inputs: library fastq, genome fasta

```shell
mbarq map -f <> -g <> -a <>
```

### Output files:

maplib_demo.barcode_map.annotated: final library map with annoations
maplib_demo.barcode_map: final library map without annotations
maplib_demo.blastn: blast output for each barcode: host sequence 
maplib_demo.fasta: fasta files of barcodes and host sequences (>barcode\nhostsequence)
maplib_demo.output.bed: bedtools intersection of gff and barcode locations
maplib_demo.temp.bed: need to clean this up after completion
tnseq2_mapping.log: log file, will only have errors in it


3. Run test suit

```
pytest tests/unit/test_mapping.py -v
```

### Run `tnseq2 count` on test data
1. Display 
```
tnseq2 count --help
```

2. Run on test data
    Required Inputs:
    Optional Inputs: mapping file 
    Output Files:
        count_demo_counts_mapped: barcodes and counts, for barcodes found in the mapping file
        count_demo_counts_unmapped: barcodes and counts for barcodes not found in the mapping file
        If not mapping file is provided, all barcode and counts will be in the count_demo_counts_mapped
        
```
tnseq2 count -f tests/test_files/dnaid2023_12_test.fasta -n count_demo -o tests/test_data -m tests/test_files/ref/library_13_1.barcode_map.annotated.csv
```

3. 

```
tnseq2 count -f tests/test_files/L -n WISH_count_demo -o tests/test_data -tn GGAGGTTCACAATGTGGGAGGTCA:40:0:after

```


### `tnseq2 merge`

 Input: directory with counts for different samples
 Output: 1 csv files with with all the sample counts to be used for the analysis
 
### `tnseq2 analyze`

Input: sample counts csv, expreimental design table, control tags 
Output: log2 FC for each gene according to experimental design, z-scores/pval relative to control tags. 


## Steps:

0. (optional) Demultiplexing
1. Mapping
2. Barcode Counting
3. Analysis
 
 
## Installation

```
conda env create -f tnseq_environment.yaml
# mamba env create -f tnseq_environment.yaml
conda activate tnseq2
pip install -e . 

```


## Quick Start



## Demultiplex

## Map


## Count

### What you need:

## Analysis

## UNDER CONSTRUCTION
