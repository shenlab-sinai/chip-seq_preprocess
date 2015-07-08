#! /bin/bash

## to run diffRepeat in one command for all uniqNmultimapped.fastq files
## $1: diffRepeats output directory path
## $2: alignment directory path where unNmultimapped.fastq files are located
## $3: repbase database required for diffRepeats
## $4: output file

CWD=`pwd`
cd ${1}

cp $2/*unNmultimapped.fastq .

FQFILES=`ls *unNmultimapped.fastq`

echo "Running: diffRepeats.pl --repbase ${3} --fq $FQFILES --tbl ${4}"

diffRepeats.pl --repbase ${3} --fq $FQFILES --tbl ${4}

rm *unNmultimapped.fastq *unNmultimapped.bam
cd ${CWD}
exit 0

