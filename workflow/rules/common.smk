from pathlib import Path
import pandas as pd


data = []
with open(config["sumstats_path"], "r") as fp:
    lines = fp.readlines()

for line in lines:
    p = Path(line.strip())
    data.append((p.stem, str(p)))

analytes = (
    pd.DataFrame.from_records(data, columns=["seqid", "sumstat_path"])
    .set_index("seqid", drop=False)
    .sort_index()
)


def get_sumstats(wildcards):
    return analytes.loc[wildcards.seqid, "sumstat_path"]


def get_final_output():
    final_output = expand(
        "results/outputs/{seqid}.regenie.tsv.gz.tbi", seqid=analytes.seqid
    )

    final_output.append("results/if/inflation_factors_table.tsv")

    final_output.extend(
        expand("results/ldsc/{seqid}/{seqid}_ldsc.log", seqid=analytes.seqid)
    )

    final_output.extend(
        expand("results/metal/{seqid}/{seqid}.metal.tsv.gz", seqid=analytes.seqid)
    )

    return final_output
