rule convert_sumstats_to_metal:
    input:
        "results/pickle/{seqid}.pkl",
    output:
        "results/metal/{seqid}/{seqid}.metal.tsv.gz",
    conda:
        "../scripts/gwaspipe/environment.yml"
    params:
        format=config.get("params").get("convert_sumstats_to_metal").get("input_format"),
        config_file=config.get("params").get("convert_sumstats_to_metal").get("config_file"),
    resources:
        runtime=30
    shell:
        "python workflow/scripts/gwaspipe/src/gwaspipe.py "
        "-f {params.format} "
        "-c {params.config_file} "
        "-i {input}"
