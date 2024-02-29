

rule munge_sumstats:
    input:
        ws_path("results/pickle/{seqid}.pkl"),
    output:
        ws_path("results/ldsc/{seqid}/{seqid}.ldsc.tsv.gz"),
    conda:
        "../scripts/gwaspipe/environment.yml"
    params:
        format=config.get("params").get("munge_sumstats").get("input_format"),
        config_file=config.get("params").get("munge_sumstats").get("config_file"),
    resources:
        runtime=lambda wc, attempt: attempt * 20,
    shell:
        "python workflow/scripts/gwaspipe/src/gwaspipe.py "
        "-f {params.format} "
        "-c {params.config_file} "
        "-i {input}"


rule compute_ldscore:
    input:
        smstat=ws_path("results/ldsc/{seqid}/{seqid}.ldsc.tsv.gz"),
    output:
        ldsc=ws_path("results/ldsc/{seqid}/{seqid}_ldsc.log"),
    conda:
        "../envs/ldsc.yaml"
    resources:
        runtime=lambda wc, attempt: attempt * 60,
    params:
        ofile=lambda wildcards, output: output.ldsc.replace(".log", ""),
        ldref=config["ldscore_reference"],
    shell:
        "ldsc.py --h2 {input.smstat} "
        "--ref-ld-chr {params.ldref}/chr@ "
        "--w-ld-chr {params.ldref}/chr@ --out {params.ofile}"
