from pathlib import Path
import yaml


# Check whether configured output formats are supported by the workflow
OUTPUT_FORMATS = config.get("output_formats", ["tsv.gz"])
SUPPORTED_OUTPUT_FORMATS = {"tsv.gz", "parquet"}
unsupported_formats = (set(OUTPUT_FORMATS) - SUPPORTED_OUTPUT_FORMATS)

if unsupported_formats:
    raise WorkflowError(
        "Unsupported output format(s): "
        + ", ".join(sorted(unsupported_formats))
    )

# Map each harmonization option to the corresponding harmonization config
HARMONIZATION_OPTIONS = {
    "pre_filtering_and_harmonization": "harmonize_sumstats_pre_filtering",
    "harmonization": "harmonize_sumstats",
    "harmonization_and_post_filtering": "harmonize_sumstats_post_filtering",
}
active_options = [harm_opt for harm_opt in HARMONIZATION_OPTIONS if config["run"].get(harm_opt, False)]

if len(active_options) > 1:
    raise WorkflowError(
        "Only one harmonization mode can be enabled at a time: "
        + ", ".join(active_options)
    )

# Validate configured output formats with harmonization config
if active_options:
    active_option = active_options[0]
    params_key = HARMONIZATION_OPTIONS[active_option]
    harmonization_config_path = Path(config["params"][params_key]["config_file"])

    with harmonization_config_path.open() as f:
        harmonization_config = yaml.safe_load(f)

    run_sequence = harmonization_config["run_sequence"]
    run_steps = [step for _, step in run_sequence]

    # TSV is always required
    if "write_tsv" not in run_steps:
        raise WorkflowError(
            f"'write_tsv' is missing from run_sequence in "
            f"{harmonization_config_path}"
        )

    if "parquet" in OUTPUT_FORMATS:
        if "write_parquet" not in run_steps:
            raise WorkflowError(
                f"'write_parquet' is missing from run_sequence in "
                f"{harmonization_config_path}, but "
                f"output_formats contains 'parquet'."
            )
    elif "write_parquet" in run_steps:
        raise WorkflowError(
            f"'write_parquet' is present in run_sequence in "
            f"{harmonization_config_path}, but "
            f"'parquet' is not present in output_formats."
        )