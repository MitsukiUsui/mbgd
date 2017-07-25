IFS=$'\n'
dir=/home/mitsuki/altorf/mbgd/analyze
cmd=${dir}/calc_identity.sh

for strain in `cat ../data/strain.lst`
do
	outlog=${dir}/log/${strain}_ide.sgeout
	errlog=${dir}/log/${strain}_ide.sgeerr
	qsub -S /bin/bash -q all.q -o ${outlog} -e ${errlog}\
		${cmd} ${strain}
done
