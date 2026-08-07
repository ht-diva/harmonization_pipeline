rule harmonize_sumstats:
    input:
        sumstats=get_sumstats,
    output:
        sumstats=temp(
            expand(ws_path(
                "temp/{{sumstat_id}}/"
                "{{sumstat_id}}.gwaslab.{output_format}"
            ),
            output_format=OUTPUT_FORMATS)
        ),
        log=temp(ws_path("temp/{sumstat_id}/{sumstat_id}.gwaslab.log")),
    conda:
        "../envs/gwaspipe.yaml"
    params:
        format=config.get("params").get("harmonize_sumstats").get("input_format"),
        config_file=config.get("params")
        .get("harmonize_sumstats_post_filtering")
        .get("config_file"),
        output_path=config.get("workspace_path"),
        sumstats_sep=config.get("sumstats_sep"),
    shell:
        "gwaspipe "
        "-f {params.format} "
        "-c {params.config_file} "
        "-s '{params.sumstats_sep}' "
        "-i {input.sumstats} "
        "-o {params.output_path}"


rule post_filtering:
    input:
        sumstats=rules.harmonize_sumstats.output.sumstats,
        log=rules.harmonize_sumstats.output.log,
    output:
        sumstats=expand(
            ws_path(
                "outputs/{{sumstat_id}}/"
                "{{sumstat_id}}.gwaslab.{output_format}"
            ),
            output_format=OUTPUT_FORMATS,
        ),
        log=ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.log"),
    conda:
        "../envs/filtering.yaml"
    params:
        snpid2filter=config.get("snpid2filter"),
        filter_snpid_col=config.get("filter_snpid_col"),
        filter_keep_flag=lambda wc: "--filter_keep" if config.get("filter_keep", False) else "",
        input_args=lambda wc, input: " ".join(f"-i {path}" for path in input.sumstats),
        output_args=lambda wc, output: " ".join(f"-o {path}" for path in output.sumstats),
    shell:
        "python workflow/scripts/filtering_by_snipid.py "
        "{params.input_args} "
        "{params.output_args} "
        "-f {params.snpid2filter} "
        "--input_snpid_column SNPID "
        "--filter_snpid_column {params.filter_snpid_col} "
        "{params.filter_keep_flag} && "
        "cp {input.log} {output.log}"


rule bgzip_tabix:
    input:
        ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.tsv.gz"),
    output:
        ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.tsv.gz.tbi"),
    conda:
        "../envs/bgzip_tabix.yaml"
    shell:
        "workflow/scripts/bgzip_tabix.sh {input} {threads}"


rule create_snp_mapping_table:
    input:
        sumstats=get_sumstat(),
        sumstats_post=get_sumstat_post_filtering(),
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
        "-o {params.output_path} && "
        "bash workflow/scripts/filter_snp_mapping.sh "
        "{input.sumstats_post} "
        "{output}"
