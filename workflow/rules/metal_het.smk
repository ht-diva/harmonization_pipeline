rule standardize_metal_heterogeneity_sumstats:
    input:
        sumstats=get_sumstats,
    output:
        ws_path("metal/{seqid}/{seqid}.metal_het.tsv.gz"),
        ws_path("min_P/{seqid}.nlargest.txt"),
        ws_path("if/{seqid}.if.txt"),
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


rule generate_if_table:
    input:
        expand(ws_path("if/{seqid}.if.txt"), seqid=analytes.seqid),
    output:
        ws_path("inflation_factors_table.tsv"),
    conda:
        "../envs/create_report_table.yaml"
    params:
        input_path=ws_path("if"),
    shell:
        "python workflow/scripts/create_report_table.py -i {params.input_path} -o {output}"


rule generate_min_pvalue_table:
    input:
        expand(ws_path("min_P/{seqid}.nlargest.txt"), seqid=analytes.seqid),
    output:
        ws_path("min_pvalue_table.tsv"),
    conda:
        "../envs/create_report_table.yaml"
    params:
        input_path=ws_path("min_P"),
    shell:
        "python workflow/scripts/create_report_table.py -i {params.input_path} -o {output}"
