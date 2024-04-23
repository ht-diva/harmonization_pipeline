rule standardize_metal_heterogeneity_sumstats:
    input:
        sumstats=get_sumstats,
    output:
        ws_path("metal/{seqid}/{seqid}.metal_het.tsv.gz"),
    conda:
        "../scripts/gwaspipe/environment.yml"
    params:
        format=config.get("params")
        .get("annotate_metal_heterogeneity_sumstats")
        .get("input_format"),
        config_file=config.get("params")
        .get("annotate_metal_heterogeneity_sumstats")
        .get("config_file"),
        output_path=config.get("workspace_path"),
    resources:
        runtime=lambda wc, attempt: attempt * 30,
    shell:
        "python workflow/scripts/gwaspipe/src/gwaspipe.py "
        "-f {params.format} "
        "-c {params.config_file} "
        "-i {input} "
        "-o {params.output_path}"
