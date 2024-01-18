# Analyze

## Identify differentially abundant genes between the control (the inoculum) and treatment conditions with `mbarq analyze`

### Input/Output Files

**Required Inputs**

- Count file produced by `mbarq merge`

| barcode   | Name | Sample1 | Sample2 | ... |
|:----------| :---: | :---: | :---: | :---: |
| ACCTGGTAG | geneA | 500 | 1000 | ... |
| ACCGGGGAA | geneA | 100 | 500 | ... |
 | CCCGGGAAA | geneB | 300 | 300 | ... |


- Sample data file (CSV) in the following format:

| sampleID | treatment | 
|:---------| :---: | 
| Sample1  | control |
| Sample2  | treatment1 | 
| ...      | ... | ... |


- Name of the column indicating treatment in the sample data should be specified using ``--treatment_column`` (for the example above, `` --treatment_column treatment``)
- Treatment level that should be used as a control/baseline should be specified using ``--baseline`` (for the example above, ``--baseline control``)

**Suggested Inputs**

- We highly recommend adding control strains (i.e. strains with barcodes inserted into fitness-neutral locations) to the barcode library. This greatly facilitates quality control and analysis of the data.
- If control strains are present in the library, the control barcodes can be specified with a control file using the ``--control_file`` option. 
  - In the simplest option, the control file will only contain the barcode sequences of the control strains (1 barcode per line). 
  - If different control strains were added at different concentrations, the concentration of each barcode can be specified in the second column. 
  - If control strains included strains of different genotypes (ex. wild type as well as negative control strains), the genotype can be specified in the 3rd column. 
  - Only wild-type strains will be used for quality control and analysis. This should be specified as `wt`, `WT`, or `wildtype`. 
  - The control file should be in CSV format, and contain NO header. 

| [Required] | [Optional] | [Optional] |
|:-----------|:----------:|:----------:|
 |ACCTGGGTT | 0.005 | wt |
| CCGGAAGGT | 0.001 | wt | 


**Output Files**

- ``mbarq_merged_counts_batch.txt``: Information on sample and batch
- ``mbarq_merged_counts.correlations.csv``: Correlation for each batch
- ``mbarq_merged_counts_rra_results.csv``: Information for each gene about number of barcodes, LFC and false discovery rate

For each comparison:
- ``mbarq_merged_counts_cond1_vs_cond0.gene_summary.txt``: Summary for each gene
- ``mbarq_merged_counts_cond1_vs_cond0.report.Rmd``: MAGeCK Comparison Report
- ``mbarq_merged_counts_cond1_vs_cond0.sgrna_summary.txt``: Summary for each sgRNA


### Example Usage

```bash 

mbarq analyze -i <count_file> -s <sample_data_file> -c <control_file> \ 
--treatment_column treatment --batch_column batch --baseline control 

```

### All Options

```
mbarq analyze
Usage: mbarq analyze <options>

Options:
  -i, --count_file FILE    CSV file produced by `mbarq merge`
  -s, --sample_data FILE   CSV file containing sample data
  -c, --control_file FILE  control barcode file, see documentation for proper
                           format
  -g, --gene_name STR      column in the count file containing gene
                           identifiers [Name]
  --treatment_column STR   column in sample data file indicating treatment
  --baseline STR           treatment level to use as control/baseline, ex.
                           day0
  -n, --name STR           experiment name, by default will try to use count
                           file name
  -o, --out_dir DIR        Output directory
  -h, --help               Show this message and exit.

```
