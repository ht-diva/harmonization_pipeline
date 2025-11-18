import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import click
from pathlib import Path

@click.command()
@click.option("-i", "--input_path", required=True, help="Input path")
@click.option("-o", "--output_path", required=True, help="Output path")
def main(input_path, output_path):

    df = pd.read_csv(
        input_path,
        sep="\t",
        compression="gzip",
        engine="pyarrow",
    )

    table = pa.Table.from_pandas(df, preserve_index=False)

    pq.write_table(
        table,
        output_path,
        compression="zstd",
        use_dictionary=True,
        coerce_timestamps="ms",
        flavor="spark",
    )



if __name__ == "__main__":
    main()
