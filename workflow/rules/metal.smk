rule convert_sumstats_to_metal:
    input:
        ws_path("pickle/{seqid}.pkl"),
    output:
        ws_path("metal/{seqid}/{seqid}.metal.tsv.gz"),
    conda:
        "../scripts/gwaspipe/environment.yml"
    params:
        format=config.get("params").get("convert_sumstats_to_metal").get("input_format"),
        config_file=config.get("params")
        .get("convert_sumstats_to_metal")
        .get("config_file"),
        output_path=config.get("workspace_path"),
    resources:
        runtime=lambda wc, attempt: attempt * 30,
    shell:
        "python workflow/scripts/gwaspipe/src/gwaspipe.py "
        "-f {params.format} "
        "-c {params.config_file} "
        "-i {input} "
        "-o {params.output_path}"
