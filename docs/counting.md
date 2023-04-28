# Count

## Quantify barcoded strains in a sample with `mbarq count`

### Input/Output files

**Required Inputs**

- Sample sequencing file. Can be FASTQ of FASTA, commpressed is accepted. 
- Transposon construct structure

**Sugested Inputs**

- Mapping file. Can be the mapping file produced by ``mbarq map``, or any ``csv`` file, where the first column is titled ``barcode`` and contains the barcode sequences.
- For example:


    | barcode | barcodeID |
    |---------|-----------|
    | ATGCATG | Barcode1  |

- By default, will try to merge barcodes with edit distance of <= 2. This can be changed with ``-e``

**Output Files**



### Example Usage

```bash 

mbarq count  -f <sample.fastq.gz> -m <library_mapping_file.csv> \ 
-n ExperimentName -tn B17N13GTGTATAAGAGACAG

```

### All Options

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

## Merge: merge count files from multiple samples with `mbarq merge`

### Input/Output Files

**Required Inputs**

- **EITHER** a comma-separated list of count files to be merged **OR** a path to the directory containing all the count files to be merged.
- Prefix for the output file.

**Suggested Inputs**
- If you would like to keep one of the annotation (attribute) columns from the mapping file, it can be specified with `-a` option. 

**Output file**
- `csv` file, where columns are barcodes, attribute if specified, and one column containing counts for each sample

### Example Usage

```bash

mbarq merge -d <directory_with_count_files> -a locus_tag -n ExperimentName -o .

```

### All Options

```
Usage: mbarq merge <options>

Options:
  -i, --input_files FILE[,FILE]  list of mBARq count files to merge (comma
                                 separated)
  -d, --count_dir DIR            merge all files in the directory
  -o, --out_dir DIR              output directory
  -n, --name STR                 output file prefix
  -a, --attribute STR            Feature attribute to keep in the merged file
                                 (ex. ID, Name, locus_tag)
  -h, --help                     Show this message and exit.

```

