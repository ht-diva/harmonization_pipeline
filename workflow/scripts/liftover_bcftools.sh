#!/bin/bash

INPUT=$1
INPUT_SEP=$2
HG37FASTA=$3
HG38FASTA=$4
CHAINFILE=$5
OUTPUT_VCF=$6
THREADS=$7


bcftools norm --threads ${THREADS} -f ${HG37FASTA} -c s -Ou "${INPUT_VCF}" -- | \
bcftools +liftover --threads ${THREADS} --no-version -Ou -- -s ${HG37FASTA} -f ${HG38FASTA} -c ${CHAINFILE} | \
bcftools sort -Oz -o ${OUTPUT_VCF}
