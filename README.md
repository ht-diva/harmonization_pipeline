# harmonization_pipeline
The pipeline is designed to harmonize summary statistics based on GWASLab.

## Requirements
* Singularity

see also [environment.yml](environment.yml) and [Makefile](Makefile)

## Getting started

* `git clone --recurse-submodules https://github.com/ht-diva/harmonization_pipeline.git`
* `cd harmonization_pipeline`

The [config/config.yaml](config/config.yaml) file is a configuration file that sets various parameters and paths for setting up the pipeline. 
Here's a breakdown of its components that you must accord to your needs:

**Rule Execution Flags**

These flags control which parts of the pipeline should be executed.

```yaml
run:
  harmonization: True
  liftoverbcf_harmonization: False
  gwascatalog: False
  summarize: True
  delivery: False
  to_parquet: False
 ```

* _harmonization_: If set to True, the harmonization rules will be executed.
* _liftoverbcf_harmonization_: If set to False, the harmonization rules will not be executed after BCFtools liftover.
* _gwascatalog_: If set to False, the summary statistics will not be downloaded by GWAS Catalog.
* _summarize_: If set to True, the summarization rules will be executed.
* _delivery_: If set to False, the delivery rules will not be executed.
* _to_parquet_: If set to True, the harmonized summary statistics will be available in parquet format, along with .tsv.gz.

**Paths**

These paths define the locations of input files and directories where intermediate and final results will be stored. Change the paths as per your requirements.

```yaml
sumstats_path: config/sumstats_paths.tsv
sumstats_url: config/urls_gwascatalog_test.tsv
sumstats_suffix: ".gwas.regenie.gz"
sumstats_sep: "\t"
dest_path: "../destination"
workspace_path: "../results"
hg37_fasta_file_path: '/public_data/liftOver/human_g1k_v37.fasta'
hg38_fasta_file_path: '/public_data/liftOver/hg38.fa'
chain_file_path: '/public_data/liftOver/hg19ToHg38.over.chain.gz'
```

* _sumstats_path_: Path to the file containing a path to each summary statistics to process.
* _sumstats_url_: Path to the file containing a url to each summary statistics to download from GWAS Catalog (only when `gwascatalog: True`).
* _sumstats_suffix_: Common suffix of the input summary statistics files.
* _sumstats_sep_: Separator of the input summary statistics files.
* _dest_path_: Destination path for the final results.
* _workspace_path_: Path where intermediate results will be stored.
* _hg37_fasta_file_path_: Path to the file with the reference genome used in the input summary statistics, hg19 (GRCh37) (only when `liftoverbcf_harmonization: True`).
* _hg38_fasta_file_path_: Path to the file with the target genome used in the liftover, hg38 (GRCh38) (only when `liftoverbcf_harmonization: True`).
* _chain_file_path_: Path to the chain file to map the genomic coordinates from hg19 (GRCh37) to hg38 (GRCh38) during the liftover (only when `liftoverbcf_harmonization: True`).

**Common parameters**

These parameters are used across different steps of the pipeline.

```yaml
input_format: &iformat "regenie"
```

* _input_format_: Defines the input format for summary statistics and aliased as &iformat. Check the Input formats section for more details.

**Parameters for Specific Steps**

```yaml
params:
  harmonize_sumstats:
    input_format: *iformat
    config_file: "config/config_harmonize_sumstats.yml"
  snp_mapping:
    input_format: *iformat
    config_file: "config/config_snp_mapping.yml"
  summarize_sumstats:
    input_format: "gwaslab"
    config_file: "config/config_summarize_sumstats.yml"
```

* _harmonize_sumstats_: Parameters for the harmonize_sumstats rule.
  * _input_format_: Uses the aliased *iformat.
  * _config_file_: Path to the configuration file for harmonization.

* _snp_mapping_: Parameters for the snp_mapping rule.
  * _input_format_: Uses the aliased *iformat.
  * _config_file_: Path to the configuration file for SNP mapping.

* _summarize_sumstats_: Parameters for the summarize_sumstats rule.
  * _input_format_: Set to "gwaslab".
  * _config_file_: Path to the configuration file for summarization.



At the bottom of each rule-based configuration files described above, there is a filename transformation section that specifies how to extract some information from the input filenames.

For example:

```yaml
# Filename transformation, e.g.: seq.3007.7.gwas.regenie.gz
filename_mask: [ True, True, True, False, False, False]
filename_sep: '.'
```

This transformation section tells the pipeline how to split the filename into parts. The `filename_mask` list specifies which parts should be retained (True) and which should be discarded (False). 
The `filename_sep` specifies the separator used to split the filename into parts.

These configuration files are essential for setting up and running the analysis pipeline, ensuring that each step is correctly configured and that paths and formats are consistent across the workflow.


### Configuration file examples

Examples of configuration files for *BELIEVE*, *Decode*, *FinnGen*, *Genes & Health*, *GTEx*, *GWAS Catalog*, *INTERVAL*, *Meta-CHRIS-INTERVAL*, *UK Biobank* and *UKB-PPP* input data are given in the folder [examples](examples).

### Submitting the workflow
To submit the workflow to the HT HPC cluster, you can use the [submit.sbatch](submit.sbatch) script with the command `sbatch submit.sbatch`. Check the script to adapt it to your specific requirements.

#### Note on job names

The job name can now be displayed as rule name (`rule_***`) followed by the summary statistics identifier (`wildcards_***`) in the "COMMENT" field of `squeue`. Use the command:

`squeue --me --format="%.10i %.9P %.8j %.50k %.8u %.2t %.6M %.6D %.6R"`

with output (example):

| JOBID   | PARTITION | NAME                   | COMMENT                                  | USER     | ST | TIME  | NODES | NODELIST |
|---------|-----------|------------------------|------------------------------------------|----------|----|-------|-------|----------|
| 199xxxx | cpuq      | 72a9f3ce-8929-...      | rule_harmonize_sumstats_wildcards_seq123 | username | R  | mm:ss | 1     | cnodexx  |
| 199xxxx | cpuq      | harmonization_pipeline | (null)                                   | username | R  | mm:ss | 1     | cnodexx  |

## Input formats

Possible input formats for summary statistics (see [formatbook.json](https://github.com/ht-diva/gwaspipe/blob/dev/data/formatbook.json) for more options to add):
* *finngen*
* *vcf*
* *decode*
* *gwaslab*
* *regenie*
* *regenie_gene*
* *fastgwa*
* *ldsc*
* *fuma*
* *pickle*
* *metal_het*
* *ssf_custom*
* *gwascatalog_hm_custom*
* *gtex*

### Configuration files

This pipeline requires 6 configuration files in the folder [config](config): the main configuration file [config/config.yaml](config/config.yaml), and 5 rule-based configuration files where to specify the parameters of each step of the rule.

Examples of configuration files for *BELIEVE*, *Decode*, *FinnGen*, *Genes & Health*, *GTEx*, *GWAS Catalog*, *INTERVAL*, *Meta-CHRIS-INTERVAL*, *UK Biobank* and *UKB-PPP* input data are given in the folder [examples](examples).

## Rules description
* **harmonize_sumstats** (`harmonization: True`): <br />
*Purpose:*  Performs GWASLab harmonization on input data without filtering.<br />
*Output*: *{sumstat_id}.gwaslab.tsv.gz*: Standardized and aligned GWAS summary statistics.<br />

* **bgzip_tabix** (`harmonization: True`): <br />
*Purpose*: Creates a region-based index (CHR and POS columns) of GWAS harmonized data for fast queries.<br />
*Output*: *{sumstat_id}.gwaslab.tsv.gz.tbi*: Index of GWAS harmonized data.<br />

* **to_parquet** (`to_parquet: True`): <br />
*Purpose:*  Convert harmonized summary statistics from *.tsv.gz* to *parquet* format.<br />
*Output*: *{sumstat_id}.gwaslab.parquet*: Standardized and aligned GWAS summary statistics in *parquet* format.<br />

* **create_snp_mapping_table** (`harmonization: True`): <br />
*Purpose:* Creates mapping table to match input data and harmonized summary statistics.<br />
*Output*: *table.snp_mapping.tsv.gz*: Table that links input SNPID (and rsID when available) to harmonized SNPID. Note that in case `gwascatalog: True`, SNPID mapping is generated for each summary statistics.<br />

* **convert_to_vcf** (`liftoverbcf_harmonization: True`): <br />
*Purpose:* Convert input data to VCF, i.e. the format required for BCFtools liftover (**liftover_bcftools** ).<br />
*Output*: *{sumstat_id}.vcf*: Summary statistics in VCF format.<br />

* **liftover_bcftools** (`liftoverbcf_harmonization: True`): <br />
*Purpose:* Performs BCFtools liftover on the VCF-converted input before performing harmonization (**harmonize_sumstats**), indexing (**bgzip_tabix**), and SNPID mapping (**create_snp_mapping_table**).<br />
*Output*: *{sumstat_id}.liftover.vcf.gz*: Liftover summary statistics in VCF format.<br />

* **gwascatalog_wget** (`gwascatalog: True`): <br />
*Purpose:* Downloads summary statistics from GWAS Catalog before performing harmonization (**harmonize_sumstats**), indexing (**bgzip_tabix**), and per-sumstat_id SNPID mapping (**create_snp_mapping_table**).<br />
*Output*: *{sumstat_id}.h.tsv.gz* and *{sumstat_id}.format.txt*: Summary statistics from GWAS Catalog and format (*gwascatalog_hm* or *gwascatalog*).<br />

* **summarize_sumstats**, **quality_check**, **create_if_table**, **create_min_pvalue_table** and **create_quality_check_table**  (`summarize: True`): <br />
*Purpose*: Creates summary reports and plots of harmonized data.<br />
*Outputs*:<br />
*{sumstat_id}.png*: Includes a Manhattan plot of -log10(p-values) by chromosome/position, and a QQ plot of observed -log10(p-values) vs. expected, with thresholds for genome-wide significance.<br />
*min_pvalue_table.tsv*: Table with top association hits (SNPs with the smallest p-value in the GWAS summary statistics).<br />
*inflation_factors_table.tsv*: Table with genomic inflation factors (lambda GC, Median and Maximum chi-squared statistics).<br />
*quality_check_table.tsv*: Table with quality check information: nr. of log errors, difference of line nr. between input and harmonized summary statistics, log messages for nr. of removed variants due to bad statistics and/or during liftover, log message of (in)consistent variants.<br />

* **sync_outputs_folder**, **sync_plots** and **sync_tables**  (`delivery: True`): <br />
*Purpose*: Copies GWAS indexes, and summary reports and plots to destination folder `dest_path`.<br />
*Outputs*: Copies of *{sumstat_id}.gwaslab.tsv.gz.tbi*, *{sumstat_id}.{sumstat_id}.png*, *min_pvalue_table.tsv*, *inflation_factors_table.tsv*, and *table.snp_mapping.tsv.gz*.<br />

* **sync_snp_mapping_folder** (`gwascatalog: True`): <br />
*Purpose*: Copies per-sumstat_id SNPID mapping to destination folder `dest_path`.<br />
*Outputs*: Copies of *{sumstat_id}.table.snp_mapping.tsv.gz*.<br />

### Standardization and Harmonization

Standardization and Harmonization includes the following steps:

* Check SNP identifiers (SNPID/rsID).
* Fix chromosome notation (CHR), basepair positions (POS) and alleles (EA and NEA).
* Sanity check on statistics.
* Infer genome reference build version.
* Align alleles sorting them alphabetically. Precisely EA is the effect allele for which BETA is estimated, and EA is alphabetically lower than NEA (in case, flip the alleles to match the alphabetical order).
* Flip allele-specific statistics for mismatches: BETA = - BETA; Z = - Z; EAF = 1 - EAF.
* Build SNPID column (CHR:POS:EA:NEA) (Optional with `fixid: True` and `overwrite: True` to specify in rule-based condiguration files, basic_check step).
* Re-name and re-order columns based on GWASLab format.

See also the [GWASLab website](https://cloufield.github.io/gwaslab/).

## DAGs
Check the dags for:
* the [default](dag_default.svg) option<br />
* the [BCFtools liftover](dag_liftoverbcftools.svg) option<br />
* the [GWAS Catalog](dag_gwascatalog.svg) option<br />
* with the [delivery](dag_delivery.svg) option<br />
