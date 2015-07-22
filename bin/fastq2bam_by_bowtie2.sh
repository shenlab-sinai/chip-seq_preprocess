#! /usr/bin/env bash
## To align fastq file to genome by bowtie
## $1: fastq file
## $2: bowtie2 index file with path
## $3: parseAln.pl mapq cutoff
## $4: parseAln.pl diff cutoff
## $5: target folder
## $6: if input are pair end data
## $7[ option ]: threads to be used in alignment, default is 4

FILE=$1
BOWTIE_INDEX=$2
MAPQ=$3
DIFF=$4
TARGET=$5
PE=$6
FILENAME=$(basename "$FILE")
FQDIR=$(dirname "$FILE")
EXT="${FILENAME##*.}"
FILENAME_BASE="${FILENAME%.*}"

SAM=${FQDIR}/${FILENAME_BASE}.sam

if [ -n "$7" ]; then
    CORES=$7
else
    CORES=4
fi

if [[ "$PE" == "no" ]]; then
	case "$EXT" in
		fq | fastq | FQ | FASTQ ) bowtie2 -p ${CORES} -x ${BOWTIE_INDEX} \
									${FILE} > ${SAM} 
		    ;;
		gz | GZ ) zcat ${FILE} | zcat ${FILE} | bowtie2 -p ${CORES} -x ${BOWTIE_INDEX} \
									- > ${SAM}
		    ;;
	esac
else
  FILE2=${FQDIR}/${FILENAME/R1/R2}
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

if [[ "$PE" == "no" ]]; then
   parseAln.pl ${SAM} $MAPQ $DIFF

   samtools merge ${SAM/sam/unNmultimapped.bam} ${SAM/sam/unmapped.bam} ${SAM/sam/multimapped.bam}
   bamToFastq -i ${SAM/sam/unNmultimapped.bam} -fq ${SAM/sam/unNmultimapped.fastq}
   
   mv ${SAM/sam/unmapped.bam} $TARGET
   mv ${SAM/sam/multimapped.bam} $TARGET
   mv ${SAM/sam/uniqmapped.bam} $TARGET
   mv ${SAM/sam/unmapped.bam.bai} $TARGET
   mv ${SAM/sam/multimapped.bam.bai} $TARGET
   mv ${SAM/sam/uniqmapped.bam.bai} $TARGET
   mv ${SAM/sam/unNmultimapped.fastq} $TARGET

   rm ${SAM/sam/unNmultimapped.bam} 
fi
#todo: implement similar steps for paired end

samtools view -Sb ${SAM} > ${SAM/sam/nonSorted.bam}
samtools sort -m 5G ${SAM/sam/nonSorted.bam} ${FQDIR}/${FILENAME_BASE}
samtools index ${SAM/sam/bam}
rm ${SAM} ${SAM/sam/nonSorted.bam}
mv ${SAM/sam/bam}* $TARGET
echo "Job Completed!!!"
