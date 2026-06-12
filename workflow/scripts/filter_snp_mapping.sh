#!/bin/bash

SUMSTATS_POST="$1"
SNP_MAP="$2"

TMP=$(mktemp --suffix=.tsv.gz)

awk -F'\t' '
NR==FNR {
    if (FNR==1) {
        for (i=1;i<=NF;i++) if ($i=="SNPID") c=i
        next
    }
    keep[$c]
    next
}
FNR==1 {
    for (i=1;i<=NF;i++) if ($i=="PREVIOUS_ID") p=i
    print
    next
}
$p in keep
' \
<(zcat -f "$SUMSTATS_POST") \
<(zcat -f "$SNP_MAP") \
| gzip > "$TMP"

mv "$TMP" "$SNP_MAP"
