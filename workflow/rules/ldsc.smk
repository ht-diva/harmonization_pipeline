

rule export_to_ldsc:
    input:
        "results/pickle/{sample}.pkl",
    output:
        "results/ldsc/{sample}.ldsc.tsv.gz",
    conda:
        "../scripts/gwaspipe/environment.yml"
    params:
        format=config.get("params").get("export_to_ldsc").get("input_format"),
        config_file=config.get("params").get("export_to_ldsc").get("config_file"),
    resources:
        runtime=60
    shell:
        "python workflow/scripts/gwaspipe/src/gwaspipe.py "
        "-f {params.format} "
        "-c {params.config_file} "
        "-i {input}"


