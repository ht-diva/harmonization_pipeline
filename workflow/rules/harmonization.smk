rule harmonize_sumstats:
    input:
        sumstats=get_sumstats,
    output:
        ws_path("pickle/{seqid}.pkl"),
        #ws_path("outputs/{seqid}/{seqid}.regenie.tsv.gz"),
        ws_path("if/report_if_{seqid}.txt"),
        ws_path("plots/{seqid}.png"),
    conda:
        "../scripts/gwaspipe/environment.yml"
    params:
        format=config.get("params").get("harmonize_sumstats").get("input_format"),
        config_file=config.get("params").get("harmonize_sumstats").get("config_file"),
        output_path=config.get("workspace_path"),
    resources:
        runtime=lambda wc, attempt: attempt * 60,
    shell:
        "python workflow/scripts/gwaspipe/src/gwaspipe.py "
        "-f {params.format} "
        "-c {params.config_file} "
        "-i {input.sumstats} "
        "-o {params.output_path}"


rule create_if_table:
    input:
        expand(ws_path("if/report_if_{seqid}.txt"), seqid=analytes.seqid),
    output:
        ws_path("inflation_factors_table.tsv"),
    conda:
        "../envs/create_inflation_factors_table.yaml"
    params:
        input_path=ws_path("if"),
    shell:
        "python workflow/scripts/create_inflation_factors_table.py -i {params.input_path} -o {output}"


rule save_min_pvalue:
    input:
        ws_path("pickle/{seqid}.pkl"),
    output:
        ws_path("min_P/{seqid}.nsmallest.tsv"),
    conda:
        "../scripts/gwaspipe/environment.yml"
    params:
        format=config.get("params").get("save_min_pvalue").get("input_format"),
        config_file=config.get("params").get("save_min_pvalue").get("config_file"),
        output_path=config.get("workspace_path"),
    resources:
        runtime=lambda wc, attempt: attempt * 5,
    shell:
        "python workflow/scripts/gwaspipe/src/gwaspipe.py "
        "-f {params.format} "
        "-c {params.config_file} "
        "-i {input} "
        "-o {params.output_path}"
