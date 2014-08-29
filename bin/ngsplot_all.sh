#! /bin/bash

## to run ngs.plot in one plots with all bam
## $1: path to bams
## $2: genome of ngs.plot
## $3: project name
## $4: cores to be used for ngs.plot
## $5: fagment length (option, default=150nt)

CWD=`pwd`
cd ${1}
if [ -f config.cov.txt ]; then rm config.cov.txt; fi
for i in `ls *.bam`; do
    echo -e "${i}\t-1\t\"${i/.bam/}\"" >> config.cov.txt
done

if [ -z "${4}" ]; then ${4}=2; fi
if [ -z "${5}" ]; then ${5}=150; fi

ngs.plot.r -G ${2} -C config.cov.txt -R genebody -O ${3}_${5} -FL ${5} -N 0.33 -P ${4}
ngs.plot.r -G ${2} -C config.cov.txt -R tss -O ${3}_tss_${5} -FL ${5} -P ${4}
ngs.plot.r -G ${2} -C config.cov.txt -R tes -O ${3}_tes_${5} -FL ${5} -P ${4}

IF_INPUT_EXIST=`ls | grep "[I|i]nput.*\.bam$" | wc -l`
zero=0

if [ $IF_INPUT_EXIST -eq $zero ]; then
	echo "There is no DNA input.";
else
	genNormedNgsplotConfig.py ${1} > config.cov.norm.txt
	ngs.plot.r -G ${2} -C config.cov.norm.txt -R genebody -O ${3}_norm_${5} -FL ${5} -N 0.33 -P ${4}
	ngs.plot.r -G ${2} -C config.cov.norm.txt -R tss -O ${3}_norm_tss_${5} -FL ${5} -P ${4}
	ngs.plot.r -G ${2} -C config.cov.norm.txt -R tes -O ${3}_norm_tes_${5} -FL ${5} -P ${4}
fi

cd ${CWD}
exit 0

