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

## Quick Start

### Map

```

mbarq map -f <library_R1.fastq.gz> -g <host.fasta> -a <host.gff> -l 100 \ 
-n LibraryName -tn B17N13GTGTATAAGAGACAG

```

### Count

```

mbarq count  -f <sample.fastq.gz> -m <library_mapping_file.csv> \ 
-n ExperimentName -tn B17N13GTGTATAAGAGACAG

```


### Merge

```

mbarq merge -d <directory_with_count_files> -a locus_tag -n ExperimentName -o .

```

### Analyze

```

mbarq analyze -i <count_file> -s <sample_data_file> -c <control_file> --treatement_column treatement \
--batch_column batch --baseline control 

```

## Documentation and Walkthoughs:

Please see [mBARq documentation]() for detailed instructions and tutorials

