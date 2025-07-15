rule summarize_sumstats:
    input:
        ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.tsv.gz"),
    output:
        temp(ws_path("if/{sumstat_id}.if.txt")),
        temp(ws_path("min_P/{sumstat_id}.nlargest.txt")),
        ws_path("plots/{sumstat_id}.png"),
    conda:
        "../envs/gwaspipe.yaml"
    params:
        format=config.get("params").get("summarize_sumstats").get("input_format"),
        config_file=config.get("params").get("summarize_sumstats").get("config_file"),
        output_path=config.get("workspace_path"),
    shell:
        "gwaspipe "
        "-f {params.format} "
        "-c {params.config_file} "
        "-i {input} "
        "-o {params.output_path}"


rule create_if_table:
    input:
        expand(ws_path("if/{sumstat_id}.if.txt"), sumstat_id=analytes.sumstat_id),
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
        expand(ws_path("min_P/{sumstat_id}.nlargest.txt"), sumstat_id=analytes.sumstat_id),
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
        sumstats_sep=config.get("sumstats_sep"),
    shell:
        "gwaspipe "
        "-f {params.format} "
        "-c {params.config_file} "
        "-s '{params.sumstats_sep}' "
        "-i {input.sumstats} "
        "-o {params.output_path}"
