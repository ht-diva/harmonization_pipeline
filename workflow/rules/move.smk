rule mv_harmonized_sumstats_to_final_destination:
    input:
        ws_path("outputs/{seqid}/{seqid}.regenie.tsv.gz"),
        ws_path("outputs/{seqid}/{seqid}.regenie.tsv.gz.tbi"),
        ws_path("plots/{seqid}.png"),
        ws_path("min_P/{seqid}/{seqid}.nsmallest.tsv"),
    output:
        dst_path("outputs/{seqid}/{seqid}.regenie.tsv.gz"),
        dst_path("outputs/{seqid}/{seqid}.regenie.tsv.gz.tbi"),
        dst_path("plots/{seqid}.png"),
        dst_path("min_P/{seqid}/{seqid}.nsmallest.tsv"),
    shell:
        "rsync -avz {input} {output}"
