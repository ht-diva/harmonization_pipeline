from pathlib import Path
import pandas as pd


data = []
with open(config["samples_path"], 'r') as fp:
    lines = fp.readlines()

for line in lines:
    p = Path(line.strip())
    data.append((p.stem, str(p)))

samples = (
    pd.DataFrame.from_records(data, columns=['sample_name', 'sample_path']).set_index("sample_name", drop=False).sort_index()
)


def get_sumstats(wildcards):
    return samples.loc[wildcards.sample, "sample_path"]

def get_final_output():
    final_output = expand("results/regenie/{sample}.regenie.tsv.gz.tbi", sample=samples.sample_name)

    return final_output