#! /bin/bash
## To align fastq file to genome by bowtie
## $1: fastq file
## $2: bowtie index file with path
## $3: target folder
## $4[ option ]: threads to be used in alignment, default is 4

FILE=$1
BOWTIE_INDEX=$2
FILENAME=$(basename "$FILE")
FQDIR=$(dirname "$FILE")
EXT="${FILENAME##*.}"
FILENAME_BASE="${FILENAME%.*}"
SAM=${FQDIR}/${FILENAME_BASE}.sam
if [ -n "$4" ]; then
    CORES=$4
else
    CORES=4
fi

case "$EXT" in
	fq | fastq | FQ | FASTQ ) bowtie -p ${CORES} -q -m 1 -v 3 --sam --best --strata ${BOWTIE_INDEX} \
	                                ${FILE} > ${SAM}
	    ;;
	gz | GZ ) zcat ${FILE} | bowtie -p ${CORES} -q -m 1 -v 3 --sam --best --strata ${BOWTIE_INDEX} \
	                                - > ${SAM}
	    ;;
esac

samtools view -Sb ${SAM} > ${SAM/sam/nonSorted.bam}
samtools sort -m 5G ${SAM/sam/nonSorted.bam} ${FQDIR}/${FILENAME_BASE}
samtools index ${SAM/sam/bam}
rm ${SAM} ${SAM/sam/nonSorted.bam}
mv ${SAM/sam/bam}* ${3}

