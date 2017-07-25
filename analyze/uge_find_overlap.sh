IFS=$'\n'
dir=/home/mitsuki/altorf/mbgd/analyze
cmd=${dir}/find_overlap.sh

for strain in `cat ../data/strain.lst`
do
	outlog=${dir}/log/${strain}_ovr.sgeout
	errlog=${dir}/log/${strain}_ovr.sgeerr
	qsub -S /bin/bash -q all.q -o ${outlog} -e ${errlog}\
		${cmd} ${strain}
done
