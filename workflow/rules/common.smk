from pathlib import Path
import pandas as pd
import datetime
import yaml

# Define input for the rules
data = []
with open(config["sumstats_path"], "r") as fp:
    lines = fp.readlines()

for line in lines:
    p = Path(line.strip())
    seqid = ".".join(p.stem.split(".")[:3])
    data.append((seqid, str(p)))

analytes = (
    pd.DataFrame.from_records(data, columns=["seqid", "sumstat_path"])
    .set_index("seqid", drop=False)
    .sort_index()
)


def get_sumstats(wildcards):
    return analytes.loc[wildcards.seqid, "sumstat_path"]


def get_sumstat():
    return analytes.iloc[0]["sumstat_path"]


def ws_path(file_path):
    return str(Path(config.get("workspace_path"), file_path))


def dst_path(file_path):
    return str(Path(config.get("dest_path"), file_path))


# Load cluster config
with open("slurm/cluster_config.yaml") as f:
    cluster_config = yaml.safe_load(f)

with open("slurm/config.yaml") as f:
    main_config = yaml.safe_load(f)

def get_resources(rule_name, attempt):
    cfg = cluster_config[rule_name]
    mem_mb = cfg["mem_base"] + attempt * cfg["mem_per_attempt"]
    runtime_minutes = attempt * cfg["runtime_per_attempt"]
    runtime_str = str(datetime.timedelta(minutes=runtime_minutes))
    return {
        "threads": main_config["cores"],
        "mem_mb": mem_mb,
        "runtime": runtime_str
    }


def get_final_output():
    final_output = []

    if config.get("run").get("harmonization"):
        final_output.extend(
            expand(
                ws_path("outputs/{seqid}/{seqid}.gwaslab.tsv.gz"),
                seqid=analytes.seqid,
            )
        )
        final_output.extend(
            expand(
                ws_path("plots/{seqid}.png"),
                seqid=analytes.seqid,
            )
        )

    if config.get("run").get("summarize"):
        final_output.append(ws_path("min_pvalue_table.tsv")),
        final_output.append(ws_path("inflation_factors_table.tsv")),
        final_output.append(ws_path("snp_mapping/table.snp_mapping.tsv.gz"))

    if config.get("run").get("delivery"):
        final_output.append(dst_path("tables_delivery.done")),
        final_output.extend(
            expand(
                dst_path("plots/{seqid}.png"),
                seqid=analytes.seqid,
            )
        ),
        final_output.extend(
            expand(
                dst_path("outputs/{seqid}/.delivery.done"),
                seqid=analytes.seqid,
            )
        )

    if config.get("run").get("annotation"):
        final_output.extend(
            expand(
                ws_path("outputs/{seqid}/{seqid}.gwaslab.tsv.gz"),
                seqid=analytes.seqid,
            )
        )
        final_output.extend(
            expand(
                ws_path("outputs/{seqid}/{seqid}.gwaslab.tsv.gz.tbi"),
                seqid=analytes.seqid,
            )
        )

    return final_output
