# Identifying Salmonella colonisation determinants in mice

This analysis is based on the BarSeq screen described in [Import of aspartate and malate by DcuABC drives H2/fumarate respiration to promote initial Salmonella gut-lumen colonization in mice](https://doi.org/10.1016%2Fj.chom.2020.04.013).
In the paper, the authors generated and sequenced a barcoded library of Salmonella mutants. This library was then used to infect LCM mice, and the fecal pellets from the infected mice were collected on days 1, 2,3, and 4 post infection. Each of the fecal samples was then sequenced to count the abundances of each of the barcoded mutants. The goal of the analysis was to identify which mutants are lost (and thus which genes are important for Salmonella pathogenesis) on different days post infection. 



## Setup

1. Make sure you have followed [Installation instructions](install.md)
2. Download and unpack [the test data](walkthrough_downloads/nguyenb.tar.gz). After running the command below, you should see a directory named `nguyenb_walkthrough`, which should contain all the data you need for this walkthrough. 

```shell
tar -xvf nguyenb.tar.gz
cd nguyenb_walkthrough
ls 
```


2. Make sure `mbarq` is installed and you have created and activated the `mbarq` environment

```shell

conda activate mbarq

```


## Mapping 

First thing in the analysis of random barcode mutagenesis experiment is figuring out the position of each barcode in the host genome. This can be accomplished with `mbarq map`.

**Files needed for this analysis**:

`library_11_1_sub_1.fq.gz` : a subsample of the library sequencing file. The reads in this file will contain barcodes + host DNA sequence, and allow us to identify the location of each barcoded insertion.

`SL1344.fna` is a genome sequence of the Salmonella strain used to construct the library. `SL1344.gff` is a matching annotation file. 

```{note}

For `mbarq` to run, you need to specify the transposon construct structure used for the experiment. Specifically, you need to specify the conserved IR motif, the length of the barcode, length of the spacer (if there is any) between barcode and IR, and their relative position to each other. Here's an example of a read from the `library_11_1_sub_1.fq.gz` file. The transposon conserved sequence (IR) is shown in **bold** and the barcode (17 nt long) is shown in `color`. The spacer sequence (13 nt long) is shown in lower case, and the host sequence in upper case. For `mbarq`, this will translate into `-tn B17N13GTGTATAAGAGACAG`

``GGGACCAAAGTACTAGA``tcagggttgagat**GTGTATAAGAGACAG**ATTGTATTCGCC[...]

```

To map barcodes to insertions run:

```shell

mbarq map -f library_11_1_sub_1.fq.gz -g SL1344.fna -a SL1344.gff \
          -tn B17N13GTGTATAAGAGACAG -n nguyenb_library_map
 
```

The final results will be saved in `nguyenb_library_map.annotated.csv`. Note that this was done on a test dataset, so we did not apply any filtering to the results. In reality, we recommend filtering out barcodes that are supported by only a few reads. The filtering threshold will vary from dataset to dataset, anywhere from 10 to 100 could be reasonable (for example, specifying `-l 10` will filter out any barcodes supported by less than 10 reads). 


## Counting

Now we are ready to analyze our samples (i.e. fecal pellets collected from different mice on different days p.i.). For each sample, you would need to run `mbarq count` command to generate a table of barcode counts.


**Files needed for this analysis**:

`library_11_1_map.annotated.csv`: a full map of the mutant library (the one we generated above was done using a subsample of the data, and would not contain all the barcodes). Contains information about chromosomal location of each of the barcodes.

`dnaid1315_124_subsample.fasta.gz`: a subsample of the sequencing file generated from one of the fecal samples. The reads in this file will contain just the barcode sequences (no host DNA), and allow us to count abundance of each of the barcodes in the sample.

To count the barcodes run:

```shell

mbarq count -f dnaid1315_124_subsample.fasta.gz -m nguyenb_library_map.annotated.csv \
          -tn B17N13GTGTATAAGAGACAG 

```

You can examine the resulting count table: `dnaid1315_124_subsample_mbarq_counts.csv`


## Merging count files

After generating the count files for each of your samples, you can merge them together into a single file using `mbarq merge`. To demonstrate this, we will be using 2 previously generated count files for samples dnaid1315_17 and dnaid1315_18

**Files needed for this analysis**:

`dnaid1315_17_mbarq_counts.csv` and `dnaid1315_18_mbarq_counts.csv` are examples of barcode count files. 

To merge count tables into a single table run

```shell

mbarq merge -i dnaid1315_17_mbarq_counts.csv,dnaid1315_18_mbarq_counts.csv -a Name -n nguyenb_counts

```
You can examine the resulting count table: `nguyenb_counts_mbarq_merged_counts.csv`


```{note}
You can also place all the count files into the same directory, and specify directory name with `-d` instead of listing all the file names as was shown above.

```


## Analysis

The goal of this experiment was to identify potential fitness factors on different days of Salmonella infection. Thus we want to compare mutant abundances on each day to the inoculum (mutant library + control strains cultured in LB), labeled as `d0`, to samples from `d1`, `d2`, `d3`, and `d4`.



**Files needed for this analysis**:

`library_11_1_mbarq_merged_counts.csv` contains counts for all the samples used in the experiment. 

`sample_data.csv` contains sample metadata, i.e. day post infection and mouse ID for each of the samples. 

`control_strains.csv` contains a list of barcoded wild-type isogenic strains used as a control in the study.

```{note}
You can read more about the `sample_data.csv` and `control_strains.csv` file formats in [Analysis section of documentation](analysis.md)

```

```shell

mkdir results
mbarq analyze -i library_11_1_mbarq_merged_counts.csv -s sample_data.csv -c control_strains.csv --treatment_column day --baseline d0 -o results 

```
`mbarq analyze` creates a folder containing `library_11_1_mbarq_merged_counts_rra_results.csv` that lists log fold changes (LFC) and false discovery rates (FDR) for each gene in the library. You can upload this file to the [mBARq App](https://asintsova-mbarq-app-home-0bmqtg.streamlit.app/) to create heatmaps, perform functional analysis and visualize the results in the context of KEGG metabolic maps.





