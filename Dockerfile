FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="e64bfd7a7e9008a3232718227cf69d40c99e257067e56da97187008e376bb59a"

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
#   source: workflow/envs/gwaspipe.yaml
#   prefix: /conda-envs/6160ee07ce3509ba51bd41a874dda23a
#   name: gwaspipe
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - python<=3.10
#     - pip==24
#     - mscorefonts
#     - pip:
#         - git+https://github.com/ht-diva/gwaspipe.git@main
RUN mkdir -p /conda-envs/6160ee07ce3509ba51bd41a874dda23a
COPY workflow/envs/gwaspipe.yaml /conda-envs/6160ee07ce3509ba51bd41a874dda23a/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/6e056d31662ab0bd2fd3fba49416042f --file /conda-envs/6e056d31662ab0bd2fd3fba49416042f/environment.yaml && \
    mamba env create --prefix /conda-envs/a160f42d06f9d24b41c5cbece52b682d --file /conda-envs/a160f42d06f9d24b41c5cbece52b682d/environment.yaml && \
    mamba env create --prefix /conda-envs/eb561a2c57a89b268cbdebae74913b82 --file /conda-envs/eb561a2c57a89b268cbdebae74913b82/environment.yaml && \
    mamba env create --prefix /conda-envs/31fc19a9498faffb09aa18f9246db95e --file /conda-envs/31fc19a9498faffb09aa18f9246db95e/environment.yaml && \
    mamba env create --prefix /conda-envs/6160ee07ce3509ba51bd41a874dda23a --file /conda-envs/6160ee07ce3509ba51bd41a874dda23a/environment.yaml && \
    mamba clean --all -y