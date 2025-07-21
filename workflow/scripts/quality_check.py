import click
import gzip
import re
from pathlib import Path

@click.command()
@click.option("-gs", "--gwas_sumstats_path", required=True, help="GWAS summary statistics path")
@click.option("-hs", "--harm_sumstats_path", required=True, help="Harmonized summary statistics path")
@click.option("-l", "--harm_log_path", required=True, help="Harmonization log path")
@click.option("-o", "--output_path", required=True, help="Output path")
def main(gwas_sumstats_path, harm_sumstats_path, harm_log_path, output_path):
    gwas_sumstats = Path(gwas_sumstats_path)
    harm_sumstats = Path(harm_sumstats_path)
    harm_log = Path(harm_log_path)

    # Count lines of input GWAS sumstats
    # Skeep metadata if VCF file
    open_fun = gzip.open if gwas_sumstats.suffix == ".gz" else open
    if gwas_sumstats.suffix in {".vcf", ".vcf.gz"}:
        nr_input = sum(1 for line in open_fun(gwas_sumstats, "rt") if not line.startswith("#"))
    else:
        nr_input = sum(1 for _ in open_fun(gwas_sumstats, "rt")) - 1

    # Count lines of harmonized sumstats
    nr_lines = 0
    lines2 = []
    with gzip.open(harm_sumstats, "rt") as fp:
        for line in fp:
            nr_lines += 1
            line = line.rstrip("\n")
            if line:
                lines2.append(line)
                if len(lines2) > 2:
                    lines2.pop(0)
    nr_harm = max(0, nr_lines - 1)

    # Check whether the file is corrupt,
    # i.e. unexpected end with last two lines of unequal length
    input_separator = "\t"
    harm_status = "CORRUPT" if len(lines2[0].split(input_separator)) != len(lines2[1].split(input_separator)) else "OK"

    # Calculate line difference
    nr_delta = nr_input - nr_harm

    # Log information
    log_error_nr = 0
    badstat_log = []
    liftover_log = []
    inconsistent_log = []

    timestamp_pat = re.compile(r'^[\d:/\s-]+-\s*')
    with harm_log.open("r") as fp:
        for line in fp:
            line = timestamp_pat.sub("", line).strip()
            if "error" in line.lower():
                log_error_nr += 1
            if "variants with bad statistics" in line.lower():
                badstat_log.append(line)
            if "removed unmapped variants" in line.lower():
                liftover_log.append(line)
            if "not consistent" in line.lower():
                inconsistent_log.append(line)
            if "variants with inconsistent values" in line.lower():
                inconsistent_log.append(line)

    badstat_log = badstat_log or ["NONE"]
    liftover_log = liftover_log or ["NONE"]
    inconsistent_log = inconsistent_log or ["NONE"]

    header = "\t".join([
        "GWAS_SumStats",
        "Harm_SumStats",
        "Harm_Status",
        "Nr_InputLines",
        "Nr_HarmLines",
        "Nr_DeltaLines",
        "Nr_ErrorsLog",
        "Log_BadStatistics",
        "Log_Liftover",
        "Log_Inconsistency"
    ])

    data = [
        str(gwas_sumstats),
        str(harm_sumstats),
        str(harm_status),
        str(nr_input),
        str(nr_harm),
        str(nr_delta),
        str(log_error_nr),
        "; ".join(badstat_log),
        "; ".join(liftover_log),
        "; ".join(inconsistent_log)
    ]

    output = Path(output_path)
    with open(output, "w") as fp:
        fp.write(header + "\n")
        fp.write("\t".join(data) + "\n")


if __name__ == "__main__":
    main()
