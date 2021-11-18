# Analysis of RB-TNSeq Data

Allows mapping of transposon library, as well as barcode counting


## Installation:

- Create and activate conda environment 
```
git clone ...
mamba env create -f tnseq_environment.yaml
conda activate tnseq2
tnseq2 --help

```
### Run `tnseq2 maplib` on test data

1. Display help information

```
tnseq2 maplib --help

```
2. Run `maplib`
    Required inputs: library fastq, genome fasta
    Optional inputs: annotation gff
    For test dataset have to set l=0, because only a few reads are analyzed
    Output files: 
        maplib_demo.barcode_map.annotated: final library map with annoations
        maplib_demo.barcode_map: final library map without annotations
        maplib_demo.blastn: blast output for each barcode: host sequence 
        maplib_demo.fasta: fasta files of barcodes and host sequences (>barcode\nhostsequence)
        maplib_demo.output.bed: bedtools intersection of gff and barcode locations
        maplib_demo.temp.bed: need to clean this up after completion
        tnseq2_mapping.log: log file, will only have errors in it

```
tnseq2 maplib -f tests/test_files/library_13_1_1.fq -r tests/test_files/library_13_1_2.fq -a tests/test_files/ref/Salmonella_genome+plasmids.gff -g tests/test_files/ref/Salmonella_genome_FQ312003.1_SL1344.fasta --name maplib_demo -o tests/test_data -l 0

```
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
