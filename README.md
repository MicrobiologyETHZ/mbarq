# Analysis of BarSeq Data

The main steps of the workflow involve:

1. Mapping each barcode to insertion location in the genome.
2. Profiling barcode abundances across samples.
3. Mutant fitness analyses.
4. Exploratory analysis using [mBARq web app](https://microbiomics.io/tools/mbarq-app/)

## Installation:

- Clone the repository and create and activate the conda environment

```
git clone https://github.com/MicrobiologyETHZ/mbarq.git
cd mbarq
#conda env create -f mbarq_environment.yaml
mamba env create -f mbarq_environment.yaml
conda activate mbarq
pip install -e .
mbarq --help

```

## Quick Start

- Map each barcode to the insertion location in the genome

```

mbarq map -f <library_R1.fastq.gz> -g <host.fasta> -a <host.gff> -l 100 \ 
-n LibraryName -tn B17N13GTGTATAAGAGACAG

```

- Profile barcode abundances for each sample

```

mbarq count  -f <sample.fastq.gz> -m <library_mapping_file.csv> \ 
-n ExperimentName -tn B17N13GTGTATAAGAGACAG

```


- Merge barcode counts from multiple samples into the final table

```

mbarq merge -d <directory_with_count_files> -a locus_tag -n ExperimentName -o .

```

- Identify enriched/depleted genes between treatments and control

```

mbarq analyze -i <count_file> -s <sample_data_file> -c <control_file> --treatment_column treatment \
--batch_column batch --baseline control 

```

## Documentation and Walkthrough:

Please see [mBARq documentation](https://mbarq.readthedocs.io/en/latest/) for detailed instructions and tutorials.

If you use mBARq, please cite:

> **mBARq: a versatile and user-friendly framework for the analysis of DNA barcodes from transposon insertion libraries, knockout mutants, and isogenic strain populations**
> 
> Anna Sintsova, Hans-Joachim Ruscheweyh, Christopher M Field, Lilith Feer, Bidong D Nguyen, Benjamin Daniel, Wolf-Dietrich Hardt, Julia A Vorholt, Shinichi Sunagawa#
> > 
> _Bioinformatics_ (2024)
> 
> doi: [10.1093/bioinformatics/btae078](https://doi.org/10.1093/bioinformatics/btae078)

## Data and code to reproduce mBARq manuscript figures

- Use the instructions below to run the jupyter notebook used to produce figures for the [mbarq manuscript]()

```bash

git clone --branch manuscript https://github.com/MicrobiologyETHZ/mbarq.git
# or git checkout manuscript
cd mbarq/manuscript
tar -xvzf manuscript_files.tar.gz
mamba env create -f manuscript_environment.yaml
conda activate mbarq_manuscript
mkdir figures
jupyter notebook

```
