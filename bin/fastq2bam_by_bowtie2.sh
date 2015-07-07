#! /usr/bin/env bash
## To align fastq file to genome by bowtie
## $1: fastq file
## $2: bowtie2 index file with path
## $3: target folder
## $4: if input are pair end data
## $5[ option ]: threads to be used in alignment, default is 4

FILE=$1
BOWTIE_INDEX=$2
FILENAME=$(basename "$FILE")
FQDIR=$(dirname "$FILE")
EXT="${FILENAME##*.}"
FILENAME_BASE="${FILENAME%.*}"
PE=$4

SAM=${FQDIR}/${FILENAME_BASE}.sam

if [ -n "$5" ]; then
    CORES=$5
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
   samtools view -Shb -f 4 ${SAM} > ${SAM/sam/unmapped.bam}
   
   samtools view -Sh -F 4 ${SAM} > ${SAM/sam/mapped.sam}
 
   samtools view -SH ${SAM} > ${SAM/sam/multimapped.sam} 
   grep "XS:i" ${SAM/sam/mapped.sam} >> ${SAM/sam/multimapped.sam}
   samtools view -Shb ${SAM/sam/multimapped.sam} > ${SAM/sam/multimapped.bam}
 
   samtools view -SH ${SAM} > ${SAM/sam/uniqmapped.sam} 
   grep "AS:i" ${SAM/sam/mapped.sam} | grep -v "XS:i" >> ${SAM/sam/uniqmapped.sam}
   cp ${SAM/sam/uniqmapped.sam} ${SAM}
   samtools view -Shb ${SAM/sam/uniqmapped.sam} > ${SAM/sam/uniqmapped.bam}

   samtools merge ${SAM/sam/uniqNmultimapped.bam} ${SAM/sam/uniqmapped.bam} ${SAM/sam/multimapped.bam}
   bamToFastq -i ${SAM/sam/uniqNmultimapped.bam} -fq ${SAM/sam/uniqNmultimapped.fastq}

   mv ${SAM/sam/unmapped.bam} ${3}
   mv ${SAM/sam/multimapped.bam} ${3}
   mv ${SAM/sam/uniqmapped.bam} ${3}
   mv ${SAM/sam/uniqNmultimapped.fastq} ${3}

   rm ${SAM/sam/multimapped.sam} ${SAM/sam/uniqmapped.sam} ${SAM/sam/mapped.sam} ${SAM/sam/uniqNmultimapped.bam} 
fi
#todo: implement similar steps for paired end

samtools view -Sb ${SAM} > ${SAM/sam/nonSorted.bam}
samtools sort -m 5G ${SAM/sam/nonSorted.bam} ${FQDIR}/${FILENAME_BASE}
samtools index ${SAM/sam/bam}
rm ${SAM} ${SAM/sam/nonSorted.bam}
mv ${SAM/sam/bam}* ${3}
