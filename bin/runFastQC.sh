#! /bin/bash

## to run FastQC
## $1: FastQC output directory path
## $2: num cores to use
## $3: fastq file name
## $4: whether paired end

FILE1=$3
PE=$4

fastqc -o $1 -t $2 $FILE1
if [[ "$PE" == "yes" ]]; then
    fastqc -o $1 -t $2 ${FILE1/R1/R2}
fi

