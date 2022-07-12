.. mbarq documentation master file, created by
   sphinx-quickstart on Thu Jun 16 11:33:27 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to mBARq's documentation!
=================================

Transposon mutagenesis is a powerful technique that allows identification of bacterial fitness factors under different environmental conditions. Recently, a number of studies have used  barcoded transposon mutant libraries to increase the throughput of the experiments. mBARq allows easy processing and analysis of barcoded mutant libraries for any transposon construct (Tn5, *mariner*, *etc*).


Workflow
^^^^^^^^

.. image:: images/mbarq_workflow.png

The main steps of the workflow involve:

1. Mapping of each barcode to insertion location in the genome.
2. Profiling barcode abundances across samples.
3. Mutant fitness analyses.
4. Exploratory analysis using `mBARq web app`_

.. _mBARq web app: https://share.streamlit.io/asintsova/mbarq_app/main/app.py

Installation
^^^^^^^^^^^^

You will need ``conda`` to install and run ``mBARq``.

.. important::

   It is highly recommended that you install `mamba` as it greatly speeds up the environment creation.

   .. code-block::

      conda install mamba -n base -c conda-forge


Option 1
--------

- Download :download:`this environment file <mbarq_environment_complete.yaml>` and run

.. code-block::

   mamba env create -f mbarq_environment_complete.yaml
   conda activate mbarq
   mbarq --help

Option 2
--------

- Clone the repository and create and activate conda environment

.. code-block::

   git clone https://github.com/MicrobiologyETHZ/mbarq.git
   cd mbarq
   mamba env create -f mbarq_environment.yaml
   conda activate mbarq
   pip install -e .
   mbarq --help


Quick Start
^^^^^^^^^^^

- Map each barcode to insertion location in the genome


.. code-block::

   mbarq map -f <library_R1.fastq.gz> -g <host.fasta> -a <host.gff> -l 100 \
   -n LibraryName -tn B17N13GTGTATAAGAGACAG


- Profile barcode abundances for each sample


.. code-block::

   mbarq count  -f <sample.fastq.gz> -m <library_mapping_file.csv> \
   -n ExperimentName -tn B17N13GTGTATAAGAGACAG



- Merge barcode counts from multiple samples into final table


.. code-block::

   mbarq merge -d <directory_with_count_files> -a locus_tag -n ExperimentName -o .


- Identify enriched/deplted genes between treatments and control

.. code-block::

   mbarq analyze -i <count_file> -s <sample_data_file> -c <control_file> --treatement_column treatement \
   --batch_column batch --baseline control


.. toctree::
   :maxdepth: 1
   :caption: User Guide:

   mapping
   counting
   analysis
   explore

.. toctree::
   :maxdepth: 1
   :caption: Walkthroughs:

   salmonella
   mariner

.. toctree::
   :maxdepth: 1
   :caption: Documentation:

   mbarq

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
