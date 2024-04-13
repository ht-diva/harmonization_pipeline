rule sync_tables:
    input:
        table_minp=ws_path("min_pvalue_table.tsv"),
        table_if=ws_path("inflation_factors_table.tsv"),
    output:
        touch(protected(dst_path("tables_delivery.done"))),
    params:
        table_minp=dst_path("min_pvalue_table.tsv"),
        table_if=dst_path("inflation_factors_table.tsv"),
    resources:
        runtime=lambda wc, attempt: attempt * 10,
    shell:
        """
        rsync -rlptoDvz {input.table_minp} {params.table_minp} && \
        rsync -rlptoDvz {input.table_if} {params.table_if}"""


rule sync_plots_parallel:
    input:
        files=expand(
            ws_path("plots/{file}"), file=glob_wildcards(ws_path("plots/{file}")).file
        ),
    output:
        touch(protected(dst_path("plots_delivery.done"))),
    params:
        batch_size=100,
        plots=directory(dst_path("plots/")),
    resources:
        runtime=lambda wc, attempt: attempt * 120,
    run:
        for batch in range(0, len(input.files), params.batch_size):
            batch_files = input.files[batch : batch + params.batch_size]
            batch_files_str = " ".join(batch_files)
            shell("rsync -rlptoDvz --progress {batch_files_str} {params.plots}")


rule sync_outputs_folder:
    output:
        touch(protected(dst_path("outputs_delivery.done"))),
    params:
        batch_size=100,
        folders=get_folders(ws_path("outputs")),
        output_folders=dst_path("outputs/"),
    resources:
        runtime=lambda wc, attempt: attempt * 480,
    run:
        for batch in range(0, len(input.folders), params.batch_size):
            batch_folders = input.folders[batch : batch + params.batch_size]
            batch_folders_str = " ".join(batch_folders)
            shell(
                "rsync -rlptoDvz --progress {batch_folders_str} {params.output_folders}"
            )
