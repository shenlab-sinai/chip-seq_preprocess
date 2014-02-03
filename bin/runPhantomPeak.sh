#! /bin/bash

## to run phantompeak
## $1: bam file (rmdup)

FOLDER=$(dirname ${1})
BAM=$(basename ${1})
RAN=$RANDOM
mkdir -p ${FOLDER}/${RAN}
run_spp_nodups.R -savp -rf -s=-500:5:500 -c=${1} -out=${1/.bam/.txt} \
    -tmpdir=${FOLDER}/${RAN}
rm -rf ${FOLDER}/${RAN}
