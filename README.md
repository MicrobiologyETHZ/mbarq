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

## 1. Mapping: identify insertion sites for RB-Seq library with `mbarq map`

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

```
Usage: mbarq map <options>

Options:
  -f, --forward FILE           input file for reads in forward orientation;
                               FASTQ formatted; gz is ok.  [required]
  -g, --genome FILE            reference genome in FASTA format  [required]
  -a, --gff FILE               annotation file in GFF format
  -n, --name STR               unique library name, by default will try to use
                               FASTQ filename
  -tn, --transposon STR         transposon construct structure, consisting of the following:
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
                               Note: relative position of barcode and conserved sequence motif matters, 
                               i.e. if conserved sequence motif comes before the barcode,
                               it should be written as GTGTATAAGAGACAGN13B17.
                                [B17N13GTGTATAAGAGACAG]
  -o, --out_dir DIR            output directory [.]
  -l, --filter_low_counts INT  filter out barcodes supported by [INT] or less
                               reads [0]
  -ft, --feat_type STR         feature type in the GFF file to be used for
                               annotation, e.g. gene, exon, CDS [gene]
  --attributes STR[,STR]       Feature attributes to extract from GFF file
                               [ID,Name,locus_tag]
  --closest_gene               for barcodes not directly overlapping a
                               feature, report the closest feature [False]
  -h, --help                   Show this message and exit.

```

### Re-annotate mapping file

Library map file `library.map.csv` can be re-annotated (for example, if you would like to use different feature type or attributes) without re-mapping using `mbarq annotate-mapped`. 

**Example Usage**
```shell

mbarq annotate-mapped -i <library.map.csv> -a <host.gff> -ft gene --attributes ID,Name,locus_tag

```

**All Options**

```
Usage: mbarq annotate-mapped <options>

Options:
  -i, --barcode_file FILE  unannotated barcode file produced by "map"
  -a, --gff FILE           annotation file in gff format
  -n, --name STR           unique library name, by default will try to use
                           FASTQ filename
  -o, --out_dir DIR        output directory [.]
  -ft, --feat_type STR     feature type in the gff file to be used for
                           annotation, e.g. gene, exon, CDS [gene]
  --attributes STR[,STR]   Feature attributes to extract from annotation file
                           [ID,Name,locus_tag]
  --closest_gene           for barcodes not directly overlapping a feature,
                           report the closest feature [False]
  -h, --help               Show this message and exit.

```


## Count: quantify barcoded strains in a sample with `mbarq count`

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

```
Usage: mbarq count <options>

Options:
  -f, --forward FILE       input file for reads in forward orientation; FASTQ
                           formatted; gz is ok.  [required]
  -m, --mapping_file FILE  Mapping file produced by `mbarq map`.
                           Alternatively, will accept any csv file, where  the
                           first column is titled "barcode", and contains the
                           barcode sequences.

                           Example:

                           barcode,barcodeID
                           AGACCAGTACATGACGGGTATCTCTCTGCCACTCCTGTAT,Tag_1
                           [optional]
  -o, --out_dir DIR        output directory [.]
  -n, --name STR           unique library name, by default will try to use
                           FASTQ filename
  -tn, --transposon STR     transposon construct structure, consisting of the following:
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
                           Note: relative position of barcode and conserved sequence motif matters, 
                           i.e. if conserved sequence motif comes before the barcode,
                           it should be written as GTGTATAAGAGACAGN13B17. For WISH data, use GGAGGTTCACAATGTGGGAGGTCAB40
                            [B17N13GTGTATAAGAGACAG]

  -e, --edit_distance INT  merge barcodes with edit distances <= [INT] [2]
  -h, --help               Show this message and exit.

```

