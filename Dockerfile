FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="c656609b7061aaac2e06e4657068c68e6e7f78eda29ceca74574347aed3b15a2"

RUN apt-get update > /dev/null && \
    apt-get install --no-install-recommends --yes \
    build-essential

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
#   source: workflow/scripts/gwaspipe/environment.yml
#   prefix: /conda-envs/ab67c3cfb8e1a5ad9d9cb7824966853e
#   name: gwaspipe
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - python<=3.10
#     - pip==23.3.2
#     - mscorefonts
#     - pip:
#         - gwaslab==3.4.38
#         - pandas==2.2
#         - pyarrow==15.0
#         - ruamel.yaml==0.18.5
#         - click==8.1.7
#         - loguru==0.7.2
#         - cloup==3.0.4
#         - mpmath==1.3.0
#         - pybedtools==0.9.1
RUN mkdir -p /conda-envs/ab67c3cfb8e1a5ad9d9cb7824966853e
COPY workflow/scripts/gwaspipe/environment.yml /conda-envs/ab67c3cfb8e1a5ad9d9cb7824966853e/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/6e056d31662ab0bd2fd3fba49416042f --file /conda-envs/6e056d31662ab0bd2fd3fba49416042f/environment.yaml && \
    mamba env create --prefix /conda-envs/a160f42d06f9d24b41c5cbece52b682d --file /conda-envs/a160f42d06f9d24b41c5cbece52b682d/environment.yaml && \
    mamba env create --prefix /conda-envs/ab67c3cfb8e1a5ad9d9cb7824966853e --file /conda-envs/ab67c3cfb8e1a5ad9d9cb7824966853e/environment.yaml && \
    mamba clean --all -y
