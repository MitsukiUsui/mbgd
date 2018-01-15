IFS=$'\n'
dir=/home/mitsuki/altorf/mbgd/analyze
cmd=${dir}/calc_identity.sh
#strainFilepath="strain.lst"
#strainFilepath="../data/ecoli/strain.lst"
strainFilepath="/home/mitsuki/altorf/mbgd/data/streptomyces/strain.lst"

for strain in `cat ${strainFilepath}`
do
	outlog=${dir}/log/${strain}_ide.sgeout
	errlog=${dir}/log/${strain}_ide.sgeerr
	qsub -S /bin/bash -q all.q -o ${outlog} -e ${errlog}\
		${cmd} ${strain}
done
