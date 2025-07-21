rule sync_tables:
    input:
        table_minp=ws_path("min_pvalue_table.tsv"),
        table_if=ws_path("inflation_factors_table.tsv"),
        table_qc=ws_path("quality_check_table.tsv"),
        table_snp_mapping=ws_path("snp_mapping/table.snp_mapping.tsv.gz"),
    output:
        touch(protected(dst_path("tables_delivery.done"))),
    conda:
        "../envs/delivery_sync.yaml"
    params:
        table_minp=dst_path("min_pvalue_table.tsv"),
        table_if=dst_path("inflation_factors_table.tsv"),
        table_qc=dst_path("quality_check_table.tsv"),
        table_snp_mapping=dst_path("table.snp_mapping.tsv.gz"),
    shell:
        """
        rsync -rlptoDvz {input.table_minp} {params.table_minp} && \
        rsync -rlptoDvzc {input.table_minp} {params.table_minp} && \
        rsync -rlptoDvz {input.table_snp_mapping} {params.table_snp_mapping} && \
        rsync -rlptoDvzc {input.table_snp_mapping} {params.table_snp_mapping} && \
        rsync -rlptoDvz {input.table_if} {params.table_if} && \
        rsync -rlptoDvzc {input.table_if} {params.table_if} && \
        rsync -rlptoDvz {input.table_qc} {params.table_qc} && \
        rsync -rlptoDvzc {input.table_qc} {params.table_qc}"""


rule sync_plots:
    input:
        ws_path("plots/{sumstat_id}.png"),
    output:
        protected(dst_path("plots/{sumstat_id}.png")),
    conda:
        "../envs/delivery_sync.yaml"
    shell:
        """
        rsync -rlptoDvz --progress {input} {output} && \
        rsync -rlptoDvzc {input} {output}"""


rule sync_outputs_folder:
    input:
        ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.tsv.gz.tbi"),
    output:
        touch(dst_path("outputs/{sumstat_id}/.delivery.done")),
    conda:
        "../envs/delivery_sync.yaml"
    params:
        folder=ws_path("outputs/{sumstat_id}/"),
        output_folders=dst_path("outputs/"),
    shell:
        """
        rsync -rlptoDvz --chmod "D755,F644" {params.folder} {params.output_folders} && \
        rsync -rlptoDvzc {params.folder} {params.output_folders}"""
