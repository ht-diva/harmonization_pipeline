rule summarize_sumstats:
    input:
        ws_path("outputs/{seqid}/{seqid}.gwaslab.tsv.gz"),
    output:
        temp(ws_path("if/{seqid}.if.txt")),
        temp(ws_path("min_P/{seqid}.nlargest.txt")),
        ws_path("plots/{seqid}.png"),
    conda:
        "../envs/gwaspipe.yaml"
    params:
        format=config.get("params").get("summarize_sumstats").get("input_format"),
        config_file=config.get("params").get("summarize_sumstats").get("config_file"),
        output_path=config.get("workspace_path"),
    shell:
        "python workflow/scripts/gwaspipe/src/gwaspipe/gwaspipe.py "
        "-f {params.format} "
        "-c {params.config_file} "
        "-i {input} "
        "-o {params.output_path}"


rule create_if_table:
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


rule create_min_pvalue_table:
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


rule create_snp_mapping_table:
    input:
        sumstats=get_sumstat(),
    output:
        ws_path("snp_mapping/table.snp_mapping.tsv.gz"),
    conda:
        "../envs/gwaspipe.yaml"
    params:
        format=config.get("params").get("snp_mapping").get("input_format"),
        config_file=config.get("params").get("snp_mapping").get("config_file"),
        output_path=config.get("workspace_path"),
    shell:
        "python workflow/scripts/gwaspipe/src/gwaspipe/gwaspipe.py "
        "-f {params.format} "
        "-c {params.config_file} "
        "-i {input.sumstats} "
        "-o {params.output_path}"
