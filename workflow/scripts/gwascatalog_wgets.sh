#!/bin/bash

SUMSTAT_ID=$1
MAPPING_FILE=$2
FORMAT=$3
LOG_FILE=$4

echo "$(date '+%Y/%m/%d %H:%M:%S') Starting download for $SUMSTAT_ID..." > "$LOG_FILE"

# Get the sumstat_id's paths
LINE=$(grep -P "^${SUMSTAT_ID}\t" "$MAPPING_FILE")
if [[ -z "$LINE" ]]; then
    echo "$(date '+%Y/%m/%d %H:%M:%S') ERROR: sumstat_id $SUMSTAT_ID not found in mapping file." | tee -a "$LOG_FILE"
    exit 1
fi

output_path=$(echo "$LINE" | cut -f2)
url=$(echo "$LINE" | cut -f3)
study_id=$(basename "$url")
output_dir=$(dirname "$output_path")

# Get correct file name (preferably with EFO term)
file_list=$(wget -qO- "${url}/harmonised/" | grep -o 'href="[^"]\+\.h\.tsv\.gz"' | sed 's/href="//' | sed 's/"$//')
efo_file=$(echo "$file_list" | grep -- '-EFO_' | head -n1)
non_efo_file=$(echo "$file_list" | grep -v -- '-EFO_' | head -n1)
if [[ -n "$efo_file" ]]; then
  target_basename="$efo_file"
elif [[ -n "$non_efo_file" ]]; then
  target_basename="$non_efo_file"
else
  target_basename=""
fi
if [[ -z "$target_basename" ]]; then
    echo "$(date '+%Y/%m/%d %H:%M:%S') ERROR: target summary statistics not found on GWASCatalog." | tee -a "$LOG_FILE"
    exit 1
fi

# Download summary statistics
echo "$(date '+%Y/%m/%d %H:%M:%S') Downloading from $url" | tee -a "$LOG_FILE"

wget \
  --no-parent \
  -r \
  -l1 \
  -nd \
  -e robots=off \
  --accept="$target_basename,$target_basename-meta.yaml,md5sum.txt" \
  -P "$output_dir" \
  "${url}/harmonised/"

# Check format
hfile=$(find "$output_dir" -maxdepth 1 -name '*.h.tsv.gz')
mfile=$(find "$output_dir" -maxdepth 1 -name '*-meta.yaml')
if zcat "$hfile" | head -n 1 | grep -q "hm_variant_id"; then
  echo "gwascatalog_hm_custom" > "$FORMAT"
else
  echo "ssf_custom" > "$FORMAT"
fi

# Check md5
filename=$(basename "$hfile")
expected_md5=$(grep -F " $filename" "$output_dir/md5sum.txt" | awk '{print $1}')
actual_md5=$(md5sum "$hfile" | awk '{print $1}')

if [[ "$expected_md5" == "$actual_md5" ]]; then
  echo "$(date '+%Y/%m/%d %H:%M:%S') MD5 match for $filename" | tee -a "$LOG_FILE"
else
  echo "$(date '+%Y/%m/%d %H:%M:%S') MD5 mismatch for $filename" | tee -a "$LOG_FILE"
  echo "Expected: $expected_md5" | tee -a "$LOG_FILE"
  echo "Actual: $actual_md5" | tee -a "$LOG_FILE"
fi

# Rename downloaded files to match output path
if [[ -f "$hfile" ]]; then
  mv "$hfile" "$output_path"
  mv "$mfile" "$output_dir/$SUMSTAT_ID.h.tsv.gz-meta.yaml"
  echo "$(date '+%Y/%m/%d %H:%M:%S') Renamed $(basename "$hfile") to $(basename "$output_path")" | tee -a "$LOG_FILE"
  echo "$(date '+%Y/%m/%d %H:%M:%S') Renamed $(basename "$mfile") to ${SUMSTAT_ID}.h.tsv.gz-meta.yaml" | tee -a "$LOG_FILE"
else
  echo "$(date '+%Y/%m/%d %H:%M:%S') ERROR: Expected downloaded file not found for renaming." | tee -a "$LOG_FILE"
  exit 1
fi

echo "$(date '+%Y/%m/%d %H:%M:%S') Done with $SUMSTAT_ID." | tee -a "$LOG_FILE"
