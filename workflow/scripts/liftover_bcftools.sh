#!/bin/bash

INPUT=$1
HG37FASTA=$2
HG38FASTA=$3
CHAINFILE=$4
OUTPUT_VCF=$5
THREADS=$6

bcftools norm --threads ${THREADS} -f ${HG37FASTA} -c s -Ou ${INPUT} -- | \
bcftools +liftover --threads ${THREADS} --no-version -Ou -- -s ${HG37FASTA} -f ${HG38FASTA} -c ${CHAINFILE} | \
bcftools sort -Oz -o ${OUTPUT_VCF}