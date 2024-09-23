FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="baeaa25c333875fb8e8bda73b8af453c2ecc262cf0834542e82d08b3a14abf44"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y build-essential libz-dev rsync && rm -rf /var/lib/apt/lists/*

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
#   source: workflow/envs/create_report_table.yaml
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
COPY workflow/envs/create_report_table.yaml /conda-envs/a160f42d06f9d24b41c5cbece52b682d/environment.yaml

# Conda environment:
#   source: workflow/envs/filter_infoscore.yaml
#   prefix: /conda-envs/31fc19a9498faffb09aa18f9246db95e
#   name: filter_infoscore
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - python=3.11
#     - pip==24.0
#     - pip:
#         - click==8.1.7
#         - pandas==2.2
#         - pyarrow==16.0
RUN mkdir -p /conda-envs/31fc19a9498faffb09aa18f9246db95e
COPY workflow/envs/filter_infoscore.yaml /conda-envs/31fc19a9498faffb09aa18f9246db95e/environment.yaml

# Conda environment:
#   source: workflow/scripts/gwaspipe/environment.yml
#   prefix: /conda-envs/1f62a09615f902ba65d1d0a9aee3d4aa
#   name: gwaspipe
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - python<=3.10
#     - pip==24
#     - mscorefonts
#     - pip:
#         - git+https://github.com/ht-diva/gwaslab.git@order_alleles
#         - pandas==2.2
#         - pyarrow==15.0
#         - ruamel.yaml==0.18.5
#         - click==8.1.7
#         - loguru==0.7.2
#         - cloup==3.0.4
#         - pybedtools==0.9.1
RUN mkdir -p /conda-envs/1f62a09615f902ba65d1d0a9aee3d4aa
COPY workflow/scripts/gwaspipe/environment.yml /conda-envs/1f62a09615f902ba65d1d0a9aee3d4aa/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/6e056d31662ab0bd2fd3fba49416042f --file /conda-envs/6e056d31662ab0bd2fd3fba49416042f/environment.yaml && \
    mamba env create --prefix /conda-envs/a160f42d06f9d24b41c5cbece52b682d --file /conda-envs/a160f42d06f9d24b41c5cbece52b682d/environment.yaml && \
    mamba env create --prefix /conda-envs/31fc19a9498faffb09aa18f9246db95e --file /conda-envs/31fc19a9498faffb09aa18f9246db95e/environment.yaml && \
    mamba env create --prefix /conda-envs/1f62a09615f902ba65d1d0a9aee3d4aa --file /conda-envs/1f62a09615f902ba65d1d0a9aee3d4aa/environment.yaml && \
    mamba clean --all -y
