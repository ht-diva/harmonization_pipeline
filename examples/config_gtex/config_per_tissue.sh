#!/bin/bash
# usage: ./config_per_tissue.sh Whole_Blood

TISSUE="$1"
INPUT_TSV="sumstats_from_gtex.tsv"
CONFIG_TEMPLATE="config_template.yaml"
CONFIG_OUTPUT="config.yaml"

if [ -z "$TISSUE" ]; then
    echo "Usage: $0 <TISSUE>"
    exit 1
fi

# Extract tissue-specific lines
OUTPUT_TSV="sumstats_from_gtex_${TISSUE}.tsv"
awk -F'/' -v tissue="$TISSUE" '$6 == tissue' "$INPUT_TSV" > "$OUTPUT_TSV"
echo "Created tissue-specific TSV: $OUTPUT_TSV"
echo "Number of files extracted: $(wc -l < "$OUTPUT_TSV")"

# Create tissue-specific destination folder
DEST_PATH="/project/cdh/gtex/results/${TISSUE}"
mkdir -p "$DEST_PATH"
echo "Created destination folder: $DEST_PATH"

# Copy template to new config
cp "$CONFIG_TEMPLATE" "$CONFIG_OUTPUT"

# Update config.yaml 
sed -e "s|^sumstats_path:.*|sumstats_path: config/$OUTPUT_TSV|" \
    -e "s|^dest_path:.*|dest_path: \"$DEST_PATH\"|" \
    -e "s|^workspace_path:.*|workspace_path: ../results/$TISSUE|" \
    "$CONFIG_OUTPUT" > "${CONFIG_OUTPUT}.tmp"
mv "${CONFIG_OUTPUT}.tmp" "$CONFIG_OUTPUT"

echo "Generated $CONFIG_OUTPUT for tissue $TISSUE"
