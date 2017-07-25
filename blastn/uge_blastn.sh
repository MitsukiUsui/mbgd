IFS=$'\n'
dir=/home/mitsuki/altorf/mbgd
cmd=${dir}/blastn/blastn.sh

for strain in `cat ${dir}/preprocess/strain.lst`
do
	for line in `tail -n +2 ${dir}/preprocess/sampled_cluster.csv`
	do
		family=`echo ${line}|cut -d, -f 1`

		outlog=${dir}/blastn/log/${strain}_${family}.sgeout
		errlog=${dir}/blastn/log/${strain}_${family}.sgeerr

		qsub -S /bin/bash -q all.q -o ${outlog} -e ${errlog}\
			${cmd} ${strain} ${family}	
	done
done

