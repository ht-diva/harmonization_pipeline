import pandas as pd
import click
from pathlib import Path

@click.command()
@click.option("-i", "--input_path", required=True, help="Input path")
@click.option("--input_separator", default='\t', help="Input separator")
@click.option("--input_compression", default=False, is_flag=True, help="Input compression")
@click.option("--input_snpid_column", default='SNP', help="Input SNPID column")
@click.option("-o", "--output_path", required=True, help="Output path")
@click.option("-f", "--filter_path", required=True, help="Snipid to filter path")
@click.option("--filter_snpid_column", default='ID', help="filter SNPID column")
def main(input_path, input_separator, input_compression, input_snpid_column, output_path, filter_path, filter_snpid_column):
    input = Path(input_path)
    snpid = Path(filter_path)

    snpid_df = pd.read_csv(snpid, sep='\t')
    if input_compression:
        sumstat_df = pd.read_csv(input, sep=input_separator, compression='gzip')
    else:
        sumstat_df = pd.read_csv(input, sep=input_separator)

    sumstat_df = sumstat_df[~sumstat_df[input_snpid_column].isin(snpid_df[filter_snpid_column])]

    sumstat_df.to_csv(output_path, sep='\t', compression='gzip', index=False)




if __name__ == "__main__":
    main()
