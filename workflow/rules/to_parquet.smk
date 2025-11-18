rule to_parquet:
    input:
        ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.tsv.gz"),
    output:
        ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.parquet"),
    conda:
        "../envs/filtering.yaml"
    shell:
        "workflow/scripts/to_parquet.py "
        "-i {input} "
        "-o {output} "
