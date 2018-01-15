IFS=$'\n'
dir=/home/mitsuki/altorf/mbgd/analyze
cmd=${dir}/find_overlap.sh
#strainFilepath="../data/ecoli/strain.lst"
strainFilepath="/home/mitsuki/altorf/mbgd/data/streptomyces/strain.lst"

for strain in `cat ${strainFilepath}`
do
	outlog=${dir}/log/${strain}_ovr.sgeout
	errlog=${dir}/log/${strain}_ovr.sgeerr
	qsub -S /bin/bash -q all.q -o ${outlog} -e ${errlog}\
		${cmd} ${strain}
done
