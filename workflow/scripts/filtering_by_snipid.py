import pandas as pd
import click
from pathlib import Path

@click.command()
@click.option("-i", "--input_path", required=True, help="Input path")
@click.option("--input_separator", default='\t', help="Input separator")
@click.option("--input_compression", default=False, is_flag=True, help="Input compression")
@click.option("-o", "--output_path", required=True, help="Output path")
@click.option("-f", "--filter_path", required=True, help="Snipid to filter path")
def main(input_path, input_separator, input_compression, output_path, filter_path):
    input = Path(input_path)
    snpid = Path(filter_path)

    snpid_df = pd.read_csv(snpid, sep='\t')
    if input_compression:
        sumstat_df = pd.read_csv(input, sep=input_separator, compression='gzip')
    else:
        sumstat_df = pd.read_csv(input, sep=input_separator)

    sumstat_df = sumstat_df[~sumstat_df['SNPID'].isin(snpid_df['SNPID'])]

    sumstat_df.to_csv(output_path, sep='\t', compression='gzip', index=False)




if __name__ == "__main__":
    main()
