FROM mambaorg/micromamba:latest

USER root
RUN apt-get update && apt-get install -y --no-install-recommends procps && \
    rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

RUN micromamba install -y -n base -c conda-forge git && \
    micromamba clean --all --yes

# Install mamba from conda-forge (now without defaults)
RUN micromamba install -y mamba -n base -c conda-forge && \
    micromamba clean -afy

# Clone mbarq
WORKDIR /home/mambauser
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=${MAMBA_ROOT_PREFIX}/bin:${PATH}
RUN git clone https://github.com/MicrobiologyETHZ/mbarq.git
WORKDIR /home/mambauser/mbarq

# Use bash so 'conda' works as expected in RUN steps
SHELL ["/bin/bash", "-c"]

# Create env and install mbarq
# RUN micromamba env create -f mbarq_environment.yaml && \
#     source activate mbarq && \
#     pip install -e . && \
#     conda clean -afy
RUN micromamba env create -f mbarq_environment.yaml && \
    micromamba run -n mbarq pip install -e . && \ 
    micromamba clean --all --yes

RUN micromamba run -n mbarq pip install 'setuptools<=70'

# Make the mbarq env default on PATH
ENV PATH=/opt/conda/envs/mbarq/bin:${PATH}

# Sanity check: fail build if mbarq is not callable
RUN mbarq --help >/dev/null


ENTRYPOINT ["/bin/bash"]
CMD []
