# pqtl_pipeline

## Requirements
see [environment.yml](environment.yml) and [Makefile](Makefile)

## Getting started

* git clone --recurse-submodules https://github.com/ht-diva/pqtl_pipeline.git
* cd pqtl_pipeline
* adapt the [submit.sbatch](submit.sbatch) and [config/config.yaml](config/config.yaml) files to your environment
* sbatch submit.sbatch

## DAG
check the [dag](dag.svg)

## Credits
The fine mapping branch of this pipeline has been realized by [Michele Filosi](https://github.com/filosi), see more at [https://github.com/EuracBiomedicalResearch/finemap_pipeline](https://github.com/EuracBiomedicalResearch/finemap_pipeline)
