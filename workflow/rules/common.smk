from pathlib import Path
import pandas as pd

# Define input for the rules
data = []
with open(config["sumstats_path"], "r") as fp:
    lines = fp.readlines()

for line in lines:
    p = Path(line.strip())
    sumstat_id = p.name.replace(config["sumstats_suffix"], "")
    data.append((sumstat_id, str(p)))

analytes = (
    pd.DataFrame.from_records(data, columns=["sumstat_id", "sumstat_path"])
    .set_index("sumstat_id", drop=False)
    .sort_index()
)


def get_sumstats(wildcards):
    return analytes.loc[wildcards.sumstat_id, "sumstat_path"]


def get_sumstat():
    return analytes.iloc[0]["sumstat_path"]


def ws_path(file_path):
    return str(Path(config.get("workspace_path"), file_path))


def dst_path(file_path):
    return str(Path(config.get("dest_path"), file_path))


def get_final_output():
    final_output = []

    if (
        config.get("run").get("harmonization")
        or config.get("run").get("pre_filtering_and_harmonization")
        or config.get("run").get("harmonization_and_post_filtering")
    ):
        final_output.extend(
            expand(
                ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.tsv.gz"),
                sumstat_id=analytes.sumstat_id,
            )
        )
        final_output.extend(
            expand(
                ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.tsv.gz.tbi"),
                sumstat_id=analytes.sumstat_id,
            )
        )

    if config.get("run").get("summarize"):
        final_output.append(ws_path("min_pvalue_table.tsv")),
        final_output.append(ws_path("inflation_factors_table.tsv")),
        final_output.append(ws_path("snp_mapping/table.snp_mapping.tsv.gz"))
        final_output.extend(
            expand(
                ws_path("plots/{sumstat_id}.png"),
                sumstat_id=analytes.sumstat_id,
            )
        )

    if config.get("run").get("delivery"):
        final_output.append(dst_path("tables_delivery.done")),
        final_output.extend(
            expand(
                dst_path("plots/{sumstat_id}.png"),
                sumstat_id=analytes.sumstat_id,
            )
        ),
        final_output.extend(
            expand(
                dst_path("outputs/{sumstat_id}/.delivery.done"),
                sumstat_id=analytes.sumstat_id,
            )
        )

    return final_output
