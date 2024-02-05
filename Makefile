TARGETS=create_env build_dag run

all:
	@echo "Try one of: ${TARGETS}"

create_env:
	mamba env update -n snakemake --file environment.yml

build_dag:
	snakemake --dag | dot -Tsvg > dag.svg

dry-run:
	snakemake --dry-run --profile slurm --snakefile workflow/Snakefile

run:
	snakemake --profile slurm --snakefile workflow/Snakefile
