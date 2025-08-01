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


rule quality_check:
    input:
        gwas_sumstats=get_sumstats,
        harm_sumstats=ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.tsv.gz"),
        harm_log=ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.log"),
    output:
        temp(ws_path("qc/{sumstat_id}.qc.txt")),
    conda:
        "../envs/create_report_table.yaml"
    params:
        sumstats_sep=config.get("sumstats_sep"),
    shell:
        "python workflow/scripts/quality_check.py "
        "-gs {input.gwas_sumstats} "
        "-hs {input.harm_sumstats} "
        "-l {input.harm_log} "
        "-o {output}"


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


rule create_quality_check_table:
    input:
        expand(ws_path("qc/{sumstat_id}.qc.txt"), sumstat_id=analytes.sumstat_id),
    output:
        ws_path("quality_check_table.tsv"),
    conda:
        "../envs/create_report_table.yaml"
    params:
        input_path=ws_path("qc"),
    shell:
        "python workflow/scripts/create_report_table.py -i {params.input_path} -o {output}"
