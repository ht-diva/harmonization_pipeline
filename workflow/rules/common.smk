from pathlib import Path
import pandas as pd
import json

# Store config variables for ease access
pvalcol = config["sumstat"]["pvalcol"]
pthr = config["sumstat"]["pthr"]
run_list = pd.read_csv(config["run_list"], header=0, sep="\t")

# Load genetic data information
genotype = config["genodata"]["name"]
with open(config["genodata"]["json"], "r") as f:
    gd = json.load(f)
gt = gd[genotype]

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


def get_pfile_from_chrom(wildcards):
    pfiles = gt["plinkfiles"]
    nfiles = gt["nfiles"]
    if nfiles == 1:
        ff = pfiles
    else:
        ff = [p for p in pfiles if p.find(f"dedup_{wildcards.chrom}_") > 0]
    return ff[0]


def get_sumstats(wildcards):
    return analytes.loc[wildcards.seqid, "sumstat_path"]


def ws_path(file_path):
    return str(Path(config.get("workspace_path"), file_path))


def dst_path(file_path):
    return str(Path(config.get("dest_path"), file_path))


def get_final_output():
    final_output = []

    if config.get("run").get("harmonization"):
        final_output.extend(
            expand(
                ws_path("pickle//{seqid}.pkl"),
                seqid=analytes.seqid,
            )
        )
        final_output.extend(
            expand(
                ws_path("plots/{seqid}.png"),
                seqid=analytes.seqid,
            )
        )
        final_output.append(ws_path("inflation_factors_table.tsv"))

    if config.get("run").get("save_min_pvalue"):
        final_output.extend(
            expand(
                ws_path("min_P/{seqid}.nsmallest.tsv"),
                seqid=analytes.seqid,
            )
        )

    if config.get("run").get("annotation"):
        final_output.extend(
            expand(
                ws_path("outputs/{seqid}/{seqid}.regenie.tsv.gz"),
                seqid=analytes.seqid,
            )
        )
        final_output.extend(
            expand(
                ws_path("outputs/{seqid}/{seqid}.regenie.tsv.gz.tbi"),
                seqid=analytes.seqid,
            )
        )

    if config.get("run").get("ldscore"):
        final_output.extend(
            expand(ws_path("ldsc/{seqid}/{seqid}_ldsc.log"), seqid=analytes.seqid)
        )

    if config.get("run").get("metal"):
        final_output.extend(
            expand(
                ws_path("metal/{seqid}/{seqid}.metal.tsv.gz"),
                seqid=analytes.seqid,
            )
        )

    if config.get("run").get("tiledb"):
        final_output.extend(
            expand(ws_path("vcf/{seqid}/{seqid}.vcf.gz.csi"), seqid=analytes.seqid)
        )

    if config.get("run").get("finemapping"):
        final_output.append(ws_path("all_phenos_summary.cs"))

    return final_output
