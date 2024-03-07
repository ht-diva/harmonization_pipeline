rule annotate_sumstats:
    input:
        ws_path("pickle/{seqid}.pkl"),
    output:
        ws_path("outputs/{seqid}/{seqid}.regenie.tsv.gz"),
    conda:
        "../scripts/gwaspipe/environment.yml"
    params:
        format=config.get("params").get("annotate_sumstats").get("input_format"),
        config_file=config.get("params").get("annotate_sumstats").get("config_file"),
        output_path=config.get("workspace_path"),
    resources:
        runtime=lambda wc, attempt: attempt * 60,
    shell:
        "python workflow/scripts/gwaspipe/src/gwaspipe.py "
        "-f {params.format} "
        "-c {params.config_file} "
        "-i {input} "
        "-o {params.output_path}"


rule bgzip_tabix:
    input:
        ws_path("outputs/{seqid}/{seqid}.regenie.tsv.gz"),
    output:
        ws_path("outputs/{seqid}/{seqid}.regenie.tsv.gz.tbi"),
    conda:
        "../envs/bgzip_tabix.yaml"
    shell:
        "workflow/scripts/bgzip_tabix.sh {input} {threads}"
