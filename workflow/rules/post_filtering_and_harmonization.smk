rule harmonize_sumstats:
    input:
        sumstats=get_sumstats,
    output:
        ws_path("temp/{sumstat_id}/{sumstat_id}.gwaslab.tsv.gz"),
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
        rules.harmonize_sumstats.output,
    output:
        ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.tsv.gz"),
        ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.log"),
    conda:
        "../envs/filtering.yaml"
    params:
        snpid2filter=config.get("snpid2filter"),
        filter_snpid_col=config.get("filter_snpid_col"),
    shell:
        "python workflow/scripts/filtering_by_snipid.py "
        "-i {input} "
        "-o {output} "
        "-f {params.snpid2filter} "
        "--input_snpid_column SNPID "
        "--filter_snpid_column {params.filter_snpid_col}"


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
