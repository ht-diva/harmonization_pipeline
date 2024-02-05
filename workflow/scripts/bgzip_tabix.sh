#!/bin/bash

INPUT=$1
DIRNAME=`dirname ${1}`
FILENAME=`basename ${1} .gz`


gunzip ${INPUT} && bgzip ${DIRNAME}/${FILENAME} && tabix -f -S 1 -s 1 -b 2 -e 2 ${INPUT}