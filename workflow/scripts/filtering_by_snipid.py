import pandas as pd
import gzip
import click
from pathlib import Path

@click.command()
@click.option("-i", "--input_path", required=True, multiple=True, help="Input path(s)")
@click.option("--input_separator", default='\t', help="Input separator")
@click.option("--input_snpid_column", default='SNPID', help="Input SNPID column")
@click.option("-o", "--output_path", required=True, multiple=True, help="Output path(s)")
@click.option("-f", "--filter_path", required=True, help="Snipid to filter path")
@click.option("--filter_snpid_column", default='SNPID', help="filter SNPID column")
@click.option("--filter_keep", default=False, is_flag=True, help="Whether to keep the SNPIDs or remove them")
def main(input_path, input_separator, input_snpid_column, output_path, filter_path, filter_snpid_column, filter_keep):

    snpid_df = pd.read_csv(Path(filter_path), sep='\t')

    for in_path, out_path in zip(input_path, output_path, strict=True):

        # Read input summary statistic file(s)
        in_path = Path(in_path)
        in_name = in_path.name.lower()
        if in_name.endswith(".parquet"):
            sumstat_df =pd.read_parquet(in_path, engine="pyarrow")
        elif in_name.endswith((".vcf.gz", ".vcf")):
            open_fun = gzip.open if in_path.name.lower().endswith(".gz") else open
            with open_fun(in_path, "rt") as stream:
                header_line = next(line_number for line_number, line in enumerate(stream) if line.startswith("#CHROM"))
            sumstat_df = pd.read_csv(in_path, sep="\t", compression="infer", skiprows=header_line)
        else:
            sumstat_df = pd.read_csv(in_path, sep=input_separator, compression="infer")

        # Filter by SNPID
        if filter_keep:
            sumstat_df = sumstat_df[sumstat_df[input_snpid_column].isin(snpid_df[filter_snpid_column])]
        else:
            sumstat_df = sumstat_df[~sumstat_df[input_snpid_column].isin(snpid_df[filter_snpid_column])]

        # Save filtered file
        out_path = Path(out_path)
        out_name = out_path.name.lower()
        if out_name.endswith(".parquet"):
            sumstat_df.to_parquet(out_path, engine="pyarrow", index=False)
        else:
            sumstat_df.to_csv(out_path, sep='\t', compression='gzip', index=False)




if __name__ == "__main__":
    main()
