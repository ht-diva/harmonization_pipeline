

rule munge_sumstats:
    input:
        "results/pickle/{seqid}.pkl",
    output:
        "results/ldsc/{seqid}/{seqid}.ldsc.tsv.gz",
    conda:
        "../scripts/gwaspipe/environment.yml"
    params:
        format=config.get("params").get("munge_sumstats").get("input_format"),
        config_file=config.get("params").get("munge_sumstats").get("config_file"),
    resources:
        runtime=30
    shell:
        "python workflow/scripts/gwaspipe/src/gwaspipe.py "
        "-f {params.format} "
        "-c {params.config_file} "
        "-i {input}"


rule compute_ldscore:
  input:
    smstat = 'results/ldsc/{seqid}/{seqid}.ldsc.tsv.gz',
  output:
    ldsc = 'results/ldsc/{seqid}/{seqid}_ldsc.log'
  conda:
    "../envs/ldsc.yaml"
  resources:
    runtime=60
  params:
    ofile = lambda wildcards, output: output.ldsc.replace(".log", ""),
    ldref = config['ldscore_reference']
  shell:
    "ldsc.py --h2 {input.smstat} "
    "--ref-ld-chr {params.ldref}/chr@ "
    "--w-ld-chr {params.ldref}/chr@ --out {params.ofile}"