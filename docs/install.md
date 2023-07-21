# Installation

You will need ``conda`` to install and run ``mBARq``.

```{important}
   
   It is highly recommended that you install `mamba` as it greatly speeds up the environment creation. You can read more on its installation [here](https://mamba.readthedocs.io/en/latest/installation.html).
   
```

## Option 1


- Download [this environment file](mbarq_environment_install.yaml) and run

```shell

   mamba env create -f mbarq_environment_install.yaml
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


