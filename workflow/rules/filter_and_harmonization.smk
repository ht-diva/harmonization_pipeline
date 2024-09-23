rule filter_infoscore:
    input:
        get_sumstats,
    output:
        ws_path("output/{seqid}/{seqid}.regenie.infoscore_filtered.gz"),
    conda:
        "../envs/filter_infoscore.yaml"
    params:
        snpid2filter=config.get("snpid2filter"),
    resources:
        runtime=lambda wc, attempt: attempt * 30,
    shell:
        "python workflow/scripts/filter_infoscore.py "
        "-i {input} "
        "-o {output} "
        "-f {params.snpid2filter}"


rule harmonize_sumstats:
    input:
        sumstats=rules.filter_infoscore.output,
    output:
        ws_path("pickle/{seqid}.pkl"),
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
