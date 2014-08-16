#! /bin/bash
# to rmdup of bam and index it.
# $1: bam file
# $2: output dir

TARGRET=${2}/$(basename ${1} .bam)_rmdup.bam
samtools rmdup -s ${1} ${TARGRET}
samtools index ${TARGRET}
# if [ -n "${3}" ]; then
#     samtools sort -m ${3} ${TARGRET} ${TARGRET/.bam/_sorted};
# else
#     samtools sort -m 2G ${TARGRET} ${TARGRET/.bam/_sorted};
# fi
# 
# samtools index ${TARGRET/.bam/_sorted}.bam
# mv ${TARGRET/.bam/_sorted}.bam ${TARGRET}
# mv ${TARGRET/.bam/_sorted}.bam.bai ${TARGRET}.bai
