# Credits to @filosi, adapted from
# https://github.com/EuracBiomedicalResearch/finemap_pipeline


# Divide sumstat by chromosome
rule sumstat_2_plink:
    message:
        "Separate summary stat by chromosome based on best hits"
    input:
        sumstat=ws_path("outputs/{seqid}/{seqid}.regenie.tsv.gz"),
    output:
        temp(ws_path("fm/{seqid}/tmp/chr{chrom}_sumstat.csv")),
    resources:
        runtime=lambda wc, attempt: attempt * 60,
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
        temp(ws_path("fm/{seqid}/tmp/phenotype.csv")),
    resources:
        runtime=lambda wc, attempt: attempt * 5,
    conda:
        "../envs/plink-pandas.yml"
    script:
        "../scripts/separate_pheno.py"


rule clumping:
    message:
        "Run clumping"
    input:
        smstat=ws_path("fm/{seqid}/tmp/chr{chrom}_sumstat.csv"),
    output:
        ws_path("fm/{seqid}/clumping/chr{chrom}.clumps"),
        ws_path("fm/{seqid}/clumping/chr{chrom}.log"),
    params:
        infile=get_pfile_from_chrom,
        ofile=lambda wildcards, output: output[0].replace(".clumps", ""),
        clump_logp1=config["clumping"]["logp1"],
        clump_logp2=config["clumping"]["logp2"],
        clump_r2=config["clumping"]["r2"],
        clump_kb=config["clumping"]["kb"],
        sampfile=config["sample_file"],
    resources:
        runtime=lambda wc, attempt: attempt * 100,
    conda:
        "../envs/plink2.yml"
    shell:
        """
    if [ -s {input} ]
    then
      plink2 --pfile {params.infile} --clump {input.smstat} \
      --clump-log10 'input-only' --clump-field {pvalcol} \
      --clump-log10-p1 {params.clump_logp1} --clump-log10-p2 {params.clump_logp2} \
      --clump-r2 {params.clump_r2} --clump-kb {params.clump_kb} \
      --clump-snp-field ID  --out {params.ofile} \
      --memory {resources.mem_mb} \
      --keep {params.sampfile}
    else
      touch {output[0]}
      touch {output[1]}
    fi
    """


rule enlarge_and_merge:
    message:
        "Merge overlapping clumps after enlarging to at least 1Mb"
    input:
        rules.clumping.output[0],
    output:
        ws_path("fm/{seqid}/ld/chr{chrom}.ld"),
    params:
        plinkfile=get_pfile_from_chrom,
        totsize=config["clumping"]["totsize"],
    conda:
        "../envs/plink-pandas.yml"
    resources:
        runtime=lambda wc, attempt: attempt * 100,
    script:
        "../scripts/run_susieR.py"


rule run_susieR:
    input:
        ldfile=rules.enlarge_and_merge.output,
        smstat=ws_path("fm/{seqid}/tmp/chr{chrom}_sumstat.csv"),
        phenofile=rules.cut_pheno.output,
    output:
        cs_smstat=ws_path("fm/{seqid}/cs/cs_chr{chrom}.cssmstat"),
        cs_report=ws_path("fm/{seqid}/cs/cs_chr{chrom}.csreport"),
        cs_rds=ws_path("fm/{seqid}/cs/cs_chr{chrom}_fit.rds"),  #,
    params:
        use_ld=config["susieR"]["use_ld"],
        chris_id=config["susieR"]["chris_id"],
        iter=config["susieR"]["iter"],
        min_abs_corr=config["susieR"]["min_abs_corr"],
    resources:
        runtime=lambda wc, attempt: attempt * 5,
    conda:
        "../envs/susier.yml"
    script:
        "../scripts/finemapping.R"


rule collect_by_pheno:
    input:
        expand(
            ws_path("fm/{seqid}/cs/cs_chr{chrom}.cssmstat"),
            chrom=lookup(query="pheno == '{seqid}'", within=run_list, cols="chrom"),
            seqid=analytes.seqid,
        ),
    output:
        ws_path("fm/{seqid}/summary.cs"),
    resources:
        runtime=lambda wc, attempt: attempt * 60,
    params:
        bypheno=False,
    script:
        "../scripts/aggregate_by_pheno.py"


# rule annotate_by_pheno:
#   input:
#     rules.collect_by_pheno.output
#   output:
#     ws_path("fm/{seqid}/summary_annot.cs")
#   params:
#     tophitsdir = config["sumstat"]["tophits_dir"]
#   resources:
#     mem_mb = 8000
#   script:
#     "../scripts/annotate_cs.py"


rule collect_all:
    input:
        branch(
            config["sumstat"]["annotate"] == True,
            then=expand(ws_path("fm/{seqid}/summary_annot.cs"), seqid=analytes.seqid),
            otherwise=expand(ws_path("fm/{seqid}/summary.cs"), seqid=analytes.seqid),
        ),
    output:
        ws_path("all_phenos_summary.cs"),
    resources:
        runtime=lambda wc, attempt: attempt * 60,
    params:
        bypheno=True,
    conda:
        "../envs/plink-pandas.yml"
    script:
        "../scripts/aggregate_by_pheno.py"
