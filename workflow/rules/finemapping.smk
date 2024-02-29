# Credits to @filosi
# https://github.com/EuracBiomedicalResearch/finemap_pipeline

# Store config variables for ease access
pvalcol = config["sumstat"]["pvalcol"]
pthr = config["sumstat"]["pthr"]


# Divide sumstat by chromosome
rule sumstat_2_plink:
    message:
        "Separate summary stat by chromosome based on best hits"
    input:
        sumstat=ws_path("outputs/{seqid}/{seqid}.regenie.tsv.gz"),
    output:
        temp(ws_path("fm/{seqid}/tmp/chr{chrom}_sumstat.csv")),
    resources:
        runtime=lambda wc, attempt: attempt * 100,
    params:
        pval_thr=pthr,
        pval_col=pvalcol,
    conda:
        "../envs/plink-pandas.yml"
    script:
        "../scripts/separate_gwas.py"


rule cut_pheno:
    message:
        "Extract the phenotype from the original GWAS"
    input:
        phenofile=config["pheno_file"],
    output:
        temp(ws_path("{pheno}/tmp/phenotype.csv")),
    conda:
        "../envs/plink-pandas.yml"
    script:
        "../scripts/separate_pheno.py"
