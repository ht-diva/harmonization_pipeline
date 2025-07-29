FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="3d1087e6cf2a2f39ef5c561609213911fdb414eaf67fc6f1d5a36d9afcacd3fb"

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
#   source: workflow/envs/delivery_sync.yaml
#   prefix: /conda-envs/20b7f0f77b859d9ac85875e0e8e2c471
#   name: delivery_sync
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - rsync
RUN mkdir -p /conda-envs/20b7f0f77b859d9ac85875e0e8e2c471
COPY workflow/envs/delivery_sync.yaml /conda-envs/20b7f0f77b859d9ac85875e0e8e2c471/environment.yaml

# Conda environment:
#   source: workflow/envs/filtering.yaml
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
COPY workflow/envs/filtering.yaml /conda-envs/31fc19a9498faffb09aa18f9246db95e/environment.yaml

# Conda environment:
#   source: workflow/envs/gwaspipe.yaml
#   prefix: /conda-envs/dfcee21f57d24a492fe8da5f23f9677b
#   name: gwaspipe
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - python<=3.10
#     - pip==24
#     - mscorefonts
#     - pip:
#         - git+https://github.com/ht-diva/gwaspipe.git@a0f349f
RUN mkdir -p /conda-envs/dfcee21f57d24a492fe8da5f23f9677b
COPY workflow/envs/gwaspipe.yaml /conda-envs/dfcee21f57d24a492fe8da5f23f9677b/environment.yaml

# Conda environment:
#   source: workflow/envs/liftover_bcftools.yaml
#   prefix: /conda-envs/bb7d3ca556579c4e816225676dfd5175
#   name: liftover_bcftools
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - bcftools=1.22
#     - htslib=1.22
#     - bcftools-liftover-plugin
RUN mkdir -p /conda-envs/bb7d3ca556579c4e816225676dfd5175
COPY workflow/envs/liftover_bcftools.yaml /conda-envs/bb7d3ca556579c4e816225676dfd5175/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/6e056d31662ab0bd2fd3fba49416042f --file /conda-envs/6e056d31662ab0bd2fd3fba49416042f/environment.yaml && \
    mamba env create --prefix /conda-envs/a160f42d06f9d24b41c5cbece52b682d --file /conda-envs/a160f42d06f9d24b41c5cbece52b682d/environment.yaml && \
    mamba env create --prefix /conda-envs/20b7f0f77b859d9ac85875e0e8e2c471 --file /conda-envs/20b7f0f77b859d9ac85875e0e8e2c471/environment.yaml && \
    mamba env create --prefix /conda-envs/dfcee21f57d24a492fe8da5f23f9677b --file /conda-envs/dfcee21f57d24a492fe8da5f23f9677b/environment.yaml && \
    mamba env create --prefix /conda-envs/bb7d3ca556579c4e816225676dfd5175 --file /conda-envs/bb7d3ca556579c4e816225676dfd5175/environment.yaml && \
    mamba clean --all -y
