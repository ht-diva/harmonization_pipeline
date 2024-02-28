rule mv_results_to_the_final_destination:
    input:
        "results/if/inflation_factors_table.tsv",
    output:
        reassemble_path("{input}", "final"),
    shell:
        "rsync -avz {input} {output}"
