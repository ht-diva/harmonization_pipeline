# harmonization_pipeline
The pipeline is designed to harmonize summary statistics based on GWASLab.

## Requirements
* Singularity

see also [environment.yml](environment.yml) and [Makefile](Makefile)

## Getting started

* `git clone --recurse-submodules https://github.com/ht-diva/harmonization_pipeline.git`
* `cd harmonization_pipeline`
* in [config/config.yaml](config/config.yaml):
    * adapt **run options**:
        - default option, `harmonization: True`, `summarize: True`
        - with pre- (`pre_filtering_and_harmonization: True`) or post- (`harmonization_and_post_filtering: True`) filtering option
        - with destination option (`delivery: True`)
    * adapt the **input path** to the **summary statistics** `sumstats_path`
    * adapt the **suffix** of the **summary statistics** filename `sumstats_suffix` (check filenames in `sumstats_path`)
    * adapt the **input path** to the **ID table** used to filter your data `snpid2filter` (used only with pre- or post-filtering)
    * adapt the **ID column name** of the **summary statistics** `input_snpid_col` (check files in `sumstats_path`; used only with pre-filtering)
    * adapt the **ID column name** of the **ID table** used to filter your data `filter_snpid_col` (check file at `snpid2filter`; used only with pre- or post-filtering)
    * adapt the **input_format** of `harmonize_sumstats` and `snp_mapping` based on your input data (listed in `sumstats_path`; see below for a list of possible input formats)
    * adapt the **output paths** (the output is written to the path defined by the `workspace_path`; if `delivery: True` the output is copied to `dest_path`)
* in the rule-based configuration files in [config](config), adapt the filename transformation with `filename_mask` to extract the seqid with "." separator. Examples:
    * for seq.3007.7.gwas.regenie.gz, the filename_mask is [True, True, True, False, False, False]
    * for finngen_R12_AB1_ACTINOMYCOSIS.gz, the filename_mask is [True, False]
* adapt the [submit.sbatch](submit.sbatch)
* `sbatch submit.sbatch`

### Note on job names

The job name can now be displayed as rule name in the "COMMENT" field of `squeue`. Use the command:

`squeue --me --format="%.18i %.9P %.8j %.25k %.8u %.2t %.10M %.6D %.20R"`

with output (example):

| JOBID   | PARTITION | NAME                   | COMMENT            | USER     | ST | TIME  | NODES | NODELIST |
|---------|-----------|------------------------|--------------------|----------|----|-------|-------|----------|
| 199xxxx | cpuq      | 72a9f3ce-8929-...      | harmonize_sumstats | username | R  | mm:ss | 1     | cnodexx  |
| 199xxxx | cpuq      | harmonization_pipeline | (null)             | username | R  | mm:ss | 1     | cnodexx  |

### Input formats

Possible input formats for summary statistics (see [formatbook.json](workflow/scripts/gwaspipe/data/formatbook.json) for more options to add):
* *finngen*
* *vcf*
* *decode*
* *gwaslab*
* *regenie*
* *fastgwa*
* *ldsc*
* *fuma*
* *pickle*
* *metal_het*

### Configuration files

This pipeline requires 6 configuration files in the folder [config](config): the main configuration file [config/config.yaml](config/config.yaml), and 5 rule-based configuration files where to specify the parameters of each step of the rule.

Examples of configuration files for *BELIEVE*, *CHRIS*, *Decode*, *FinnGen*, *INTERVAL* and *UKBiobank* input data are given in the folder [examples](examples).

## Rules description
* **pre_filtering** and **harmonize_sumstats** (`pre_filtering_and_harmonization: True`): <br />
*Purpose:* Filters input data (column name provided `input_snpid_col`) by an ID (SNPID or rsID) list (provided in `snpid2filter` with column name `filter_snpid_col`) and performs GWASLab harmonization on filtered data.<br />
*Output*: *{seqid}.gwaslab.tsv.gz*: Pre-filtered, standardized and aligned GWAS summary statistics.<br />

* **harmonize_sumstats** (`harmonization: True`): <br />
*Purpose:*  Performs GWASLab harmonization on input data without filtering.<br />
*Output*: *{seqid}.gwaslab.tsv.gz*: Standardized and aligned GWAS summary statistics.<br />

* **harmonize_sumstats** and **post_filtering:** (`harmonization_and_post_filtering: True`): <br />
*Purpose:* Performs GWASLab harmonization on input data and filters harmonized data by a SNPID list (provided in `snpid2filter` with column name `filter_snpid_col`).<br />
*Output*: *{seqid}.gwaslab.tsv.gz*: Standardized, aligned and post-filtered GWAS summary statistics.<br />

* **bgzip_tabix** (included in all harmonization options): <br />
*Purpose*: Creates a region-based index (CHR and POS columns) of GWAS harmonized data for fast queries.<br />
*Output*: *{seqid}.gwaslab.tsv.gz.tbi*: Index of GWAS harmonized data.<br />

* **summarize_sumstats**, **create_if_table**, **create_min_pvalue_table** and **create_snp_mapping_table**  (`summarize: True`): <br />
*Purpose*: Creates summary reports and plots of harmonized data.<br />
*Outputs*:<br />
*{seqid}.png*: Includes a Manhattan plot of -log10(p-values) by chromosome/position, and a QQ plot of observed -log10(p-values) vs. expected, with thresholds for genome-wide significance.<br />
*min_pvalue_table.tsv*: Table with top association hits (SNPs with the smallest p-value in the GWAS summary statistics).<br />
*inflation_factors_table.tsv*: Table with genomic inflation factors (lambda GC, Median and Maximum chi-squared statistics).<br />
*table.snp_mapping.tsv.gz*: Mapping file that links input SNPID (and rsID when available) to harmonized SNPID.<br />

* **sync_outputs_folder**, **sync_plots** and **sync_tables**  (`delivery: True`): <br />
*Purpose*: Copies GWAS indexes, and summary reports and plots to destination folder `dest_path`.<br />
*Outputs*: Copies of *{seqid}.gwaslab.tsv.gz.tbi*, *{seqid}.{seqid}.png*, *min_pvalue_table.tsv*, *inflation_factors_table.tsv*, and *table.snp_mapping.tsv.gz*.<br />

### GWASLab Harmonization

GWASLab Harmonization includes the following steps:

* Check SNP identifiers (SNPID/rsID).
* Fix chromosome notation (CHR), basepair positions (POS) and alleles (EA and NEA).
* Sanity check on statistics.
* Infer genome reference build version.
* Align alleles to the reference genome to ensure that alleles match the reference strand and direction (in case, flip the alleles to match the reference).
* Flip allele-specific statistics for mismatches: BETA = - BETA; Z = - Z; EAF = 1 - EAF.
* Build SNPID column (CHR:POS:NEA:EA) (Optional with `fixid: True` and `overwrite: True` to specify in rule-based condiguration files, basic_check step).
* Re-name and re-order columns based on GWASLab format.

See also the [GWASLab website](https://cloufield.github.io/gwaslab/).

## DAGs
Check the dags for:
* the [default](dag_default.svg) option<br />
* with the [pre-filtering](dag_prefiltering.svg) option, or<br />
* with the [post-filtering](dag_postfiltering.svg) option, or<br />
* with the [delivery](dag_delivery.svg) option<br />
