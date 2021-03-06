# Installation

You will need ``conda`` to install and run ``mBARq``.

```{important}
   
   It is highly recommended that you install `mamba` as it greatly speeds up the environment creation.
   
    
    conda install mamba -n base -c conda-forge
    
   
```

## Option 1


- Download [this environment file](mbarq_environment_complete.yaml) and run

```shell

   mamba env create -f mbarq_environment_complete.yaml
   conda activate mbarq
   mbarq --help
```

## Option 2

- Clone the repository and create and activate conda environment

```shell

   git clone https://github.com/MicrobiologyETHZ/mbarq.git
   cd mbarq
   mamba env create -f mbarq_environment.yaml
   conda activate mbarq
   pip install -e .
   mbarq --help

```


