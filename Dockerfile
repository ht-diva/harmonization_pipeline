FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="88db0b1103f13c7ac73a5b0a6f1ec8650dcb4159cdf95fbe3c11b4d2bd29e34c"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y build-essential libz-dev && rm -rf /var/lib/apt/lists/*

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/bgzip_tabix.yaml
#   prefix: /conda-envs/6e056d31662ab0bd2fd3fba49416042f
#   name: bgzip_tabix
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - htslib=1.16
RUN mkdir -p /conda-envs/6e056d31662ab0bd2fd3fba49416042f
COPY workflow/envs/bgzip_tabix.yaml /conda-envs/6e056d31662ab0bd2fd3fba49416042f/environment.yaml

# Conda environment:
#   source: workflow/envs/create_inflation_factors_table.yaml
#   prefix: /conda-envs/a160f42d06f9d24b41c5cbece52b682d
#   name: create_inflation_factors_table
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - python<=3.10
#     - pip==23.3.2
#     - pip:
#         - click==8.1.7
RUN mkdir -p /conda-envs/a160f42d06f9d24b41c5cbece52b682d
COPY workflow/envs/create_inflation_factors_table.yaml /conda-envs/a160f42d06f9d24b41c5cbece52b682d/environment.yaml

# Conda environment:
#   source: workflow/envs/ldsc.yaml
#   prefix: /conda-envs/9fdf8a7c84ff129c8967eb5848ece5ab
#   name: ldsc
#   channels:
#   - conda-forge
#   - bioconda
#   dependencies:
#   - python=2.7
#   - bitarray=0.8
#   - nose=1.3
#   - pybedtools=0.7
#   - ldsc==1.0.1
#   - pip
#   - pip:
#     - scipy==0.18
#     - pandas==0.20
#     - numpy==1.16
RUN mkdir -p /conda-envs/9fdf8a7c84ff129c8967eb5848ece5ab
COPY workflow/envs/ldsc.yaml /conda-envs/9fdf8a7c84ff129c8967eb5848ece5ab/environment.yaml

# Conda environment:
#   source: workflow/envs/plink-pandas.yml
#   prefix: /conda-envs/a9b8ccc53333def92898879edc6df0ca
#   name: plink-pandas
#   dependencies:
#     - bioconda::plink=1.90*
#     - bioconda::tabix=1.11
#     - python=3.11.*
#     - numpy
#     - scipy
#     - pandas
#     - pip
RUN mkdir -p /conda-envs/a9b8ccc53333def92898879edc6df0ca
COPY workflow/envs/plink-pandas.yml /conda-envs/a9b8ccc53333def92898879edc6df0ca/environment.yaml

# Conda environment:
#   source: workflow/envs/plink2.yml
#   prefix: /conda-envs/fd0f8d28d055274e6ce4cf006d07ec05
#   name: plink2
#   channels:
#     - defaults
#   dependencies:
#     - bioconda::plink2==2.00a5
RUN mkdir -p /conda-envs/fd0f8d28d055274e6ce4cf006d07ec05
COPY workflow/envs/plink2.yml /conda-envs/fd0f8d28d055274e6ce4cf006d07ec05/environment.yaml

# Conda environment:
#   source: workflow/envs/susier.yml
#   prefix: /conda-envs/be316e9967fe3846694829a2f8f63147
#   name: susier
#   channels:
#     - r
#     - conda-forge
#     - defaults
#   dependencies:
#     - r-base==4.3.0
#     - r-susier==0.12.35
#     - r-tidyverse==2.0.0
#     - r-dplyr==1.1.2
#     - r-stringr==1.5.0
#     - r-purrr==1.0.1
#     - r-glue==1.6.2
#     - r-ggplot2==3.4.2
#     - r-ggpubr==0.6.0
#     - r-data.table==1.14.8
#     - r-rcpp==1.0.*
#     - r-rcppeigen==0.3.3.*
#     - r-markdown==1.10
#     - r-rlog==0.1.0
#     - numpy==1.24.3
#     - scipy==1.9.3
#     - pandas==1.5.3
#     - pip
RUN mkdir -p /conda-envs/be316e9967fe3846694829a2f8f63147
COPY workflow/envs/susier.yml /conda-envs/be316e9967fe3846694829a2f8f63147/environment.yaml

# Conda environment:
#   source: workflow/scripts/gwaspipe/environment.yml
#   prefix: /conda-envs/d3af165ce1ac53aefab55f833f2ce3d8
#   name: gwaspipe
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - python<=3.10
#     - pip==24
#     - mscorefonts
#     - pip:
#         - gwaslab==3.4.43
#         - pandas==2.2
#         - pyarrow==15.0
#         - ruamel.yaml==0.18.5
#         - click==8.1.7
#         - loguru==0.7.2
#         - cloup==3.0.4
#         - pybedtools==0.9.1
RUN mkdir -p /conda-envs/d3af165ce1ac53aefab55f833f2ce3d8
COPY workflow/scripts/gwaspipe/environment.yml /conda-envs/d3af165ce1ac53aefab55f833f2ce3d8/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/6e056d31662ab0bd2fd3fba49416042f --file /conda-envs/6e056d31662ab0bd2fd3fba49416042f/environment.yaml && \
    mamba env create --prefix /conda-envs/a160f42d06f9d24b41c5cbece52b682d --file /conda-envs/a160f42d06f9d24b41c5cbece52b682d/environment.yaml && \
    mamba env create --prefix /conda-envs/9fdf8a7c84ff129c8967eb5848ece5ab --file /conda-envs/9fdf8a7c84ff129c8967eb5848ece5ab/environment.yaml && \
    mamba env create --prefix /conda-envs/a9b8ccc53333def92898879edc6df0ca --file /conda-envs/a9b8ccc53333def92898879edc6df0ca/environment.yaml && \
    mamba env create --prefix /conda-envs/fd0f8d28d055274e6ce4cf006d07ec05 --file /conda-envs/fd0f8d28d055274e6ce4cf006d07ec05/environment.yaml && \
    mamba env create --prefix /conda-envs/be316e9967fe3846694829a2f8f63147 --file /conda-envs/be316e9967fe3846694829a2f8f63147/environment.yaml && \
    mamba env create --prefix /conda-envs/d3af165ce1ac53aefab55f833f2ce3d8 --file /conda-envs/d3af165ce1ac53aefab55f833f2ce3d8/environment.yaml && \
    mamba clean --all -y
