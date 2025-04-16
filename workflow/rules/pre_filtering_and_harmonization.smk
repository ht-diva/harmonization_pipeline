rule pre_filtering:
    input:
        sumstats=get_sumstats,
    output:
        ws_path("temp/{seqid}/{seqid}.gwaslab.tsv.gz"),
    conda:
        "../envs/filtering.yaml"
    params:
        snpid2filter=config.get("snpid2filter"),
        input_snpid_col=config.get("input_snpid_col"),
    shell:
        "python workflow/scripts/filtering_by_snipid.py "
        "-i {input} "
        "-o {output} "
        "-f {params.snpid2filter}"
        "--input_snpid_column {params.input_snpid_col}"


rule harmonize_sumstats:
    input:
        sumstats=rules.pre_filtering.output,
    output:
        ws_path("outputs/{seqid}/{seqid}.gwaslab.tsv.gz"),
    conda:
        "../scripts/gwaspipe/environment.yml"
    params:
        format=config.get("params").get("harmonize_sumstats").get("input_format"),
        config_file=config.get("params").get("harmonize_sumstats").get("config_file"),
        output_path=config.get("workspace_path"),
    shell:
        "python workflow/scripts/gwaspipe/src/gwaspipe.py "
        "-f {params.format} "
        "-c {params.config_file} "
        "-i {input.sumstats} "
        "-o {params.output_path}"


rule bgzip_tabix:
    input:
        ws_path("outputs/{seqid}/{seqid}.gwaslab.tsv.gz"),
    output:
        ws_path("outputs/{seqid}/{seqid}.gwaslab.tsv.gz.tbi"),
    conda:
        "../envs/bgzip_tabix.yaml"
    shell:
        "workflow/scripts/bgzip_tabix.sh {input} {threads}"
