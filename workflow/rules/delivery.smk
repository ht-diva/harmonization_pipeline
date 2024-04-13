rule move_to_destination:
    input:
        table_minp=ws_path("min_pvalue_table.tsv"),
        table_if=ws_path("inflation_factors_table.tsv"),
    output:
        table_minp=dst_path("min_pvalue_table.tsv"),
        table_if=dst_path("inflation_factors_table.tsv"),
    resources:
        runtime=lambda wc, attempt: attempt * 120,
    shell:
        """
        rsync -rlptoDvz {input.table_minp} {output.table_minp} && \
        rsync -rlptoDvz {input.table_if} {output.table_if}"""
