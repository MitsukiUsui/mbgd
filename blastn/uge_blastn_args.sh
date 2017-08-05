#!/bin/bash

dir=/home/mitsuki/altorf/mbgd
cmd=${dir}/blastn/blastn.sh

sstrain=${1}
family=${2}
qstrain=${3}
outlog=${dir}/blastn/log/${sstrain}_${family}.sgeout
errlog=${dir}/blastn/log/${sstrain}_${family}.sgeerr

qsub -S /bin/bash -q all.q -o ${outlog} -e ${errlog}\
	${cmd} ${sstrain} ${family} ${qstrain}	
