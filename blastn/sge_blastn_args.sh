#!/bin/bash

dir=/home/mitsuki/altorf/mbgd
cmd=${dir}/blastn/blastn.sh
strain=${1}
family=${2}
outlog=${dir}/blastn/log/${strain}_${family}.sgeout
errlog=${dir}/blastn/log/${strain}_${family}.sgeerr

qsub -S /bin/bash -q all.q -o ${outlog} -e ${errlog}\
	${cmd} ${strain} ${family}	
