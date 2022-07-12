# Analysis of BarSeq Data

The main steps of the workflow involve:

1. Mapping each barcode to insertion location in the genome.
2. Profiling barcode abundances across samples.
3. Mutant fitness analyses.
4. Exploratory analysis using [mBARq web app](https://share.streamlit.io/asintsova/mbarq_app/main/Home.py)

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

## Quick Start

- Map each barcode to insertion location in the genome

```

mbarq map -f <library_R1.fastq.gz> -g <host.fasta> -a <host.gff> -l 100 \ 
-n LibraryName -tn B17N13GTGTATAAGAGACAG

```

- Profile barcode abundances for each sample

```

mbarq count  -f <sample.fastq.gz> -m <library_mapping_file.csv> \ 
-n ExperimentName -tn B17N13GTGTATAAGAGACAG

```


- Merge barcode counts from multiple samples into final table

```

mbarq merge -d <directory_with_count_files> -a locus_tag -n ExperimentName -o .

```

- Identify enriched/deplted genes between treatments and control

```

mbarq analyze -i <count_file> -s <sample_data_file> -c <control_file> --treatement_column treatement \
--batch_column batch --baseline control 

```

## Documentation and Walkthoughs:

Please see [mBARq documentation](https://mbarq.readthedocs.io/en/latest/) for detailed instructions and tutorials

