TARGETS=dependencies dag run unlock

all:
	@echo "Try one of: ${TARGETS}"

dependencies:
	mamba env update -n snakemake --file environment.yml

dag:
	snakemake --dag | dot -Tsvg > dag.svg

dry-run:
	snakemake --dry-run --profile slurm --snakefile workflow/Snakefile

run:
	snakemake --profile slurm --snakefile workflow/Snakefile

unlock:
	snakemake --unlock
