.. mbarq documentation master file, created by
   sphinx-quickstart on Thu Jun 16 11:33:27 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to mBARq's documentation!
=================================

.. warning::

   This documentation is currently under construction


Installation
^^^^^^^^^^^^

You will need `conda` to install and run `mBARq`.

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

mBARq Web App
^^^^^^^^^^^^^

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
