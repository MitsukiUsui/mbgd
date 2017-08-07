IFS=$'\n'
dir=/home/mitsuki/altorf/mbgd/analyze
cmd=${dir}/phase_overlap.sh
strainFilepath="../data/ecoli/strain.lst"
#strainFilepath="strain.lst"

for strain in `cat ${strainFilepath}`
do
	outlog=${dir}/log/${strain}_phs.sgeout
	errlog=${dir}/log/${strain}_phs.sgeerr
	qsub -S /bin/bash -q all.q -o ${outlog} -e ${errlog}\
		${cmd} ${strain}
done
