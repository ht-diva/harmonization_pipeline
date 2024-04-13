rule move_to_destination:
    input:
        table_minp=ws_path("min_pvalue_table.tsv"),
        table_if=ws_path("inflation_factors_table.tsv"),
        plots=ws_path("plots"),
    output:
        touch(ws_path("delivery.done")),
        table_minp=dst_path("min_pvalue_table.tsv"),
        table_if=dst_path("inflation_factors_table.tsv"),
        plots=dst_path("plots"),
    resources:
        runtime=lambda wc, attempt: attempt * 120,
    shell:
        """
        rsync -rlptoDvz {input.plots} {output.plots} && \
        rsync -rlptoDvz {input.table_minp} {output.table_minp} && \
        rsync -rlptoDvz {input.table_if} {output.table_if}"""
