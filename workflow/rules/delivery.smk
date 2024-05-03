rule sync_tables:
    input:
        table_minp=ws_path("min_pvalue_table.tsv"),
        table_if=ws_path("inflation_factors_table.tsv"),
        table_snp_mapping=ws_path("snp_mapping/table.snp_mapping.tsv.gz"),
    output:
        touch(protected(dst_path("tables_delivery.done"))),
    params:
        table_minp=dst_path("min_pvalue_table.tsv"),
        table_if=dst_path("inflation_factors_table.tsv"),
        table_snp_mapping=dst_path("table.snp_mapping.tsv.gz"),
    resources:
        runtime=lambda wc, attempt: attempt * 10,
    shell:
        """
        rsync -rlptoDvz {input.table_minp} {params.table_minp} && \
        rsync -rlptoDvz {input.table_snp_mapping} {params.table_snp_mapping} && \
        rsync -rlptoDvz {input.table_if} {params.table_if}"""


rule sync_plots:
    input:
        ws_path("plots/{seqid}.png"),
    output:
        protected(dst_path("plots/{seqid}.png")),
    resources:
        runtime=lambda wc, attempt: attempt * 30,
    shell:
        """
        rsync -rlptoDvz --progress {input} {output}"""


rule sync_outputs_folder:
    input:
        ws_path("outputs/{seqid}/{seqid}.gwaslab.tsv.gz.tbi"),
    output:
        touch(dst_path("outputs/{seqid}/.delivery.done")),
    params:
        folder=ws_path("outputs/{seqid}/"),
        output_folders=dst_path("outputs/"),
    resources:
        runtime=lambda wc, attempt: attempt * 60,
    shell:
        """
        rsync -rlptoDvz --chmod "D755,F644" {params.folder} {params.output_folders}"""
