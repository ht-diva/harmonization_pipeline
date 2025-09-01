rule gwascatalog_wget:
    input:
        sumstats_path=config.get("sumstats_path"),
    output:
        sumstats=ws_path(temp("temp/gwascatalog/{sumstat_id}/{sumstat_id}.h.tsv.gz")),
        format=ws_path(temp("temp/gwascatalog/{sumstat_id}/{sumstat_id}.format.txt")),
        log=ws_path(temp("temp/gwascatalog/{sumstat_id}/{sumstat_id}.gwascatalog_wget.log")),
    conda:
        "../envs/gwascatalog_wget.yaml"
    shell:
        "workflow/scripts/gwascatalog_wgets.sh {wildcards.sumstat_id} {input.sumstats_path} {output.format} {output.log}"


rule harmonize_sumstats:
    input:
        sumstats=get_sumstats,
        format=ws_path("temp/gwascatalog/{sumstat_id}/{sumstat_id}.format.txt"),
    output:
        ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.tsv.gz"),
        ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.log"),
    conda:
        "../envs/gwaspipe.yaml"
    params:
        config_file=config.get("params").get("harmonize_sumstats").get("config_file"),
        output_path=config.get("workspace_path"),
        sumstats_sep=config.get("sumstats_sep"),
    shell:
        "gwaspipe "
        "-f $(cat {input.format}) "
        "-c {params.config_file} "
        "-s '{params.sumstats_sep}' "
        "-i {input.sumstats} "
        "-o {params.output_path}"


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
        sumstats=get_sumstats,
        format=ws_path("temp/gwascatalog/{sumstat_id}/{sumstat_id}.format.txt"),
    output:
        ws_path("snp_mapping/{sumstat_id}/table.snp_mapping.tsv.gz"),
    conda:
        "../envs/gwaspipe.yaml"
    params:
        config_file=config.get("params").get("snp_mapping").get("config_file"),
        output_path=config.get("workspace_path"),
        sumstats_sep=config.get("sumstats_sep"),
    shell:
        "gwaspipe "
        "-f $(cat {input.format}) "
        "-c {params.config_file} "
        "-s '{params.sumstats_sep}' "
        "-i {input.sumstats} "
        "-o {params.output_path}"
