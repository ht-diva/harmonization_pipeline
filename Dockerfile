FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="234ec22e3be232aeda946a3b6a7ec095a5e60ae66f8296cb84fd32cab55912a9"

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
#   prefix: /conda-envs/eb561a2c57a89b268cbdebae74913b82
#   name: delivery_sync
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - rsync
RUN mkdir -p /conda-envs/eb561a2c57a89b268cbdebae74913b82
COPY workflow/envs/delivery_sync.yaml /conda-envs/eb561a2c57a89b268cbdebae74913b82/environment.yaml

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
#   source: workflow/scripts/gwaspipe/environment.yml
#   prefix: /conda-envs/76b872a666188b687e7041ea7992c8d8
#   name: gwaspipe
#   channels:
#     - conda-forge
#     # We want to have a reproducible setup, so we don't want default channels,
#     # which may be different for different users. All required channels should
#     # be listed explicitly here.
#     - nodefaults
#   dependencies:
#     - python=3.10.*  # or don't specify the version and use the latest stable Python
#     - pip  # pip must be mentioned explicitly, or conda-lock will fail
#     - poetry=2.*  # or 1.1.*, or no version at all -- as you want
RUN mkdir -p /conda-envs/76b872a666188b687e7041ea7992c8d8
COPY workflow/scripts/gwaspipe/environment.yml /conda-envs/76b872a666188b687e7041ea7992c8d8/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/6e056d31662ab0bd2fd3fba49416042f --file /conda-envs/6e056d31662ab0bd2fd3fba49416042f/environment.yaml && \
    mamba env create --prefix /conda-envs/a160f42d06f9d24b41c5cbece52b682d --file /conda-envs/a160f42d06f9d24b41c5cbece52b682d/environment.yaml && \
    mamba env create --prefix /conda-envs/eb561a2c57a89b268cbdebae74913b82 --file /conda-envs/eb561a2c57a89b268cbdebae74913b82/environment.yaml && \
    mamba env create --prefix /conda-envs/31fc19a9498faffb09aa18f9246db95e --file /conda-envs/31fc19a9498faffb09aa18f9246db95e/environment.yaml && \
    mamba env create --prefix /conda-envs/76b872a666188b687e7041ea7992c8d8 --file /conda-envs/76b872a666188b687e7041ea7992c8d8/environment.yaml && \
    mamba clean --all -y
