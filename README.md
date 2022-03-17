# Analysis of BarSeq Data

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

## Usage

### Identify insertion sites for RB-Seq library

**Required inputs**: 
- FASTQ file generated from sequencing randomly barcoded mutant library
- genome FASTA file of the bacteria used to generate the library
- Transposon construct structure **add diagram**.

**Suggested inputs**:
- Annotation file in GFF3 format (this will allow mapping insertion sites to genomic features). 
- Filtering parameter (``-l``, in our hands, filtering barcodes supported by less than 100 reads produced reliable library annotations. This is likely to be dataset dependent, and should be tested for each use case).
- Report closest gene (``-c``). If ``gff`` files is provided, by default, ``mbarq`` will only report features overlaping the insertion site. In addition, ``mbarq`` can report the location and distance of the closest downstream feature for barcodes that do not directly overlap any features. 

**Example Usage**

```shell

mbarq map -f <library_R1.fastq.gz> -g <host.fasta> -a <host.gff> -l 100

```

**Output files**

``library.annotated.csv``: final library map with annotations 

``library.map.csv``: final library map without annotations 

``library_mapping.log``: log file 

``library.blastn``: blast output for each barcode: host sequence;  
``library.fasta``: fasta files of barcodes and host sequences (>barcode\nhostsequence); 
``library.output.bed``: bedtools intersection of gff and barcode locations

**All Options**

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

### Quantify barcoded strains in a sample

**Required Inputs**
- Sample sequencing file. Can be FASTQ of FASTA, commpressed is accepted. 
- Transposon construct structure

**Sugested Inputs**

- Mapping file. Can be the mapping file produced by ``mbarq map``, or any ``csv`` file, where the first column is titled ``barcode`` and contains the barcode sequences. For example:

    | barcode | barcodeID |
    |---------|-----------|
    | ATGCATG | Barcode1  |

- By default, will try to merge barcodes with edit distance of <= 2. This can be changed with ``-e``

**Example Usage**
```
mbarq count  -f <sample.fastq.gz> -m <library_mapping_file.csv>
```

**All Options**

```shell

mbarq count --help
Usage: mbarq count <options>

Options:
  -f, --forward FILE       input file for reads in forward orientation; FASTQ
                           formatted; gz is ok.  [required]
  -m, --mapping_file FILE  Barcode map/annotation file in csv format.First
                           column must be titled "barcode", and contain the
                           barcode sequences. [optional]
                           
                           Example:
                           
                           barcode,barcodeID
                           AGACCAGTACATGACGGGTATCTCTCTGCCACTCCTGTAT,Tag_1
                           
  -o, --out_dir DIR        output directory [.]
  -n, --name STR           unique library name, by default will try to use
                           FASTQ filename
  -tn, --transposon STR    transposon construct structure, consisting of the following:
                              1. conserved transposon sequence, eg. GTGTATAAGAGACAG
                              2. barcode length, written as B[# of nt], eg. B17
                              3. if there are extra nucleotides between barcode and 
                                 transposon sequence, indicate with N[# of nt], eg. N13
                           Note: relative position of barcode and transpson matters, 
                           the default (RBSeq tn) represents the following construct:
                           ---|BARCODE (17 nt)|--spacer (13 nt)--|GTGTATAAGAGACAG|---HOST--
                            [B17N13GTGTATAAGAGACAG]
                            For WISH use GGAGGTTCACAATGTGGGAGGTCAB40
                            
  -e, --edit_distance INT  merge barcodes with edit distances <= [INT] [2]
  -h, --help               Show this message and exit.

```

