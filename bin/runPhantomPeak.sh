#! /bin/bash

## to run phantompeak
## $1: bam file (rmdup)
## $2: num cores to use
## $3: output directory path

FOLDER=$(dirname ${1})
BAM=$(basename ${1})
RAN=$RANDOM
OUTFILE=${3}/${BAM/.bam/.txt}
mkdir -p ${FOLDER}/${RAN}
#run_spp_nodups.R -savp -rf -s=-500:5:500 -c=${1} -out=${1/.bam/.txt} \
run_spp_nodups.R -savp -rf -s=-500:5:500 -c=${1} -out=$OUTFILE -odir=$3\
    -tmpdir=${FOLDER}/${RAN}
rm -rf ${FOLDER}/${RAN}
