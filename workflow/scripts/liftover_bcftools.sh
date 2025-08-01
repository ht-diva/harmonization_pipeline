#!/bin/bash

INPUT=$1
HG37FASTA=$2
HG38FASTA=$3
CHAINFILE=$4
OUTPUT_VCF=$5
THREADS=$6
LOG=$7

echo "$(date '+%Y/%m/%d %H:%M:%S')  - BCFtools liftover - START..."  >> "$LOG"

# Count variants before BCFtools liftover
LINES_PRE=$(bcftools view -H "$INPUT" | wc -l)
echo "$(date '+%Y/%m/%d %H:%M:%S')  - Nr. variants before BCFtools liftover: $LINES_PRE" >> "$LOG"

# BCFtools Liftover
bcftools norm --threads ${THREADS} -f ${HG37FASTA} -c s -Ou ${INPUT} -- | \
bcftools +liftover --threads ${THREADS} --no-version -Ou -- -s ${HG37FASTA} -f ${HG38FASTA} -c ${CHAINFILE} | \
bcftools sort -Oz -o "$OUTPUT_VCF"
bcftools index "$OUTPUT_VCF"

# Filter valid chromosome notations
VALID_CHR="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
bcftools view -r "$VALID_CHR" "$OUTPUT_VCF" -Oz -o "${OUTPUT_VCF}.tmp"
mv "${OUTPUT_VCF}.tmp" "$OUTPUT_VCF"
bcftools index -f "$OUTPUT_VCF"

# Count variants after BCFtools liftover
LINES_POST=$(bcftools view -H "$OUTPUT_VCF" | wc -l)
echo "$(date '+%Y/%m/%d %H:%M:%S')  - Nr. variants after BCFtools liftover: $LINES_POST" >> "$LOG"
echo "$(date '+%Y/%m/%d %H:%M:%S')  - Dropped variants during BCFtools liftover: $((LINES_PRE - LINES_POST))" >> "$LOG"

echo "$(date '+%Y/%m/%d %H:%M:%S')  - BCFtools liftover - END."  >> "$LOG"