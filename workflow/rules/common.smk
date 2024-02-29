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


def ws_path(file_path):
    return str(Path(config.get("workspace_path"), file_path))


def dst_path(file_path):
    return str(Path(config.get("dest_path"), file_path))


def get_final_output():
    final_output = []

    final_output.extend(
        expand(
            dst_path("results/outputs/{seqid}/{seqid}.regenie.tsv.gz.tbi"),
            seqid=analytes.seqid,
        )
    )

    final_output.extend(
        expand(
            dst_path("results/outputs/{seqid}/{seqid}.regenie.tsv.gz"),
            seqid=analytes.seqid,
        )
    )

    final_output.append(dst_path("results/if/inflation_factors_table.tsv"))

    if config.get("run").get("ldscore"):
        final_output.extend(
            expand(
                ws_path("results/ldsc/{seqid}/{seqid}_ldsc.log"), seqid=analytes.seqid
            )
        )

    if config.get("run").get("metal"):
        final_output.extend(
            expand(
                ws_path("results/metal/{seqid}/{seqid}.metal.tsv.gz"),
                seqid=analytes.seqid,
            )
        )

    if config.get("run").get("tiledb"):
        final_output.extend(
            expand(
                ws_path("results/vcf/{seqid}/{seqid}.vcf.gz.csi"), seqid=analytes.seqid
            )
        )

    if config.get("run").get("save_min_pvalue"):
        final_output.extend(
            expand(
                dst_path("results/min_P/{seqid}/{seqid}.nsmallest.tsv"),
                seqid=analytes.seqid,
            )
        )

    return final_output
