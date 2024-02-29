rule mv_harmonized_sumstats_to_final_destination:
    input:
        ws_path("results/outputs/{seqid}/{seqid}.regenie.tsv.gz"),
        ws_path("results/outputs/{seqid}/{seqid}.regenie.tsv.gz.tbi"),
        ws_path("results/plots/{seqid}.png"),
        ws_path("results/min_P/{seqid}/{seqid}.nsmallest.tsv"),
    output:
        dst_path("results/outputs/{seqid}/{seqid}.regenie.tsv.gz"),
        dst_path("results/outputs/{seqid}/{seqid}.regenie.tsv.gz.tbi"),
        dst_path("results/plots/{seqid}.png"),
        dst_path("results/min_P/{seqid}/{seqid}.nsmallest.tsv"),
    shell:
        "rsync -avz {input} {output}"
