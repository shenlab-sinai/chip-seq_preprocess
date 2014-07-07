#! /bin/bash
## To align fastq file to genome by bowtie
## $1: fastq file
## $2: bowtie2 index file with path
## $3: target folder
## $4[ option ]: threads to be used in alignment, default is 4

FILE=$1
BOWTIE_INDEX=$2
FILENAME=$(basename "$FILE")
FQDIR=$(dirname "$FILE")
EXT="${FILENAME##*.}"
FILENAME_BASE="${FILENAME%.*}"

## test if there is pair end fq file
FILE2=${FQDIR}/${FILENAME/R1/R2}
if [ -f ${FILE2} ]; then
	PE="Y"
	echo "${FILE2}"
else
	PE="N"
fi

SAM=${FQDIR}/${FILENAME_BASE}.sam

if [ -n "$4" ]; then
    CORES=$4
else
    CORES=4
fi

if [ "$PE" == "N" ]; then
	case "$EXT" in
		fq | fastq | FQ | FASTQ ) bowtie2 -p ${CORES} -x ${BOWTIE_INDEX} \
									${FILE} > ${SAM}
		    ;;
		gz | GZ ) zcat ${FILE} | zcat ${FILE} | bowtie2 -p ${CORES} -x ${BOWTIE_INDEX} \
									- > ${SAM}
		    ;;
	esac
else
	case "$EXT" in
		fq | fastq | FQ | FASTQ ) bowtie2 -p ${CORES} -x ${BOWTIE_INDEX} \
		                                -1 ${FILE} -2 ${FILE2} > ${SAM}
		    ;;
		gz | GZ ) mkfifo ${FILE}.fifo ${FILE2}.fifo
			zcat ${FILE} > ${FILE}.fifo & zcat ${FILE2} > ${FILE2}.fifo & \
								bowtie2 -p ${CORES} -x ${BOWTIE_INDEX} \
		            			-1 ${FILE}.fifo -2 ${FILE2}.fifo > ${SAM}
		    rm ${FILE}.fifo ${FILE2}.fifo
		    ;;
	esac
fi

samtools view -Sb ${SAM} > ${SAM/sam/nonSorted.bam}
samtools sort -m 5G ${SAM/sam/nonSorted.bam} ${FQDIR}/${FILENAME_BASE}
samtools index ${SAM/sam/bam}
rm ${SAM} ${SAM/sam/nonSorted.bam}
mv ${SAM/sam/bam}* ${3}
