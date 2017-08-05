IFS=$'\n'
dir=/home/mitsuki/altorf/mbgd
cmd=${dir}/blastn/blastn.sh

for strain in `cat ${dir}/data/ecoli/strain.lst`
do
    outlog=${dir}/blastn/log/${strain}.sgeout
    errlog=${dir}/blastn/log/${strain}.sgeerr
    qsub -S /bin/bash -q all.q -o ${outlog} -e ${errlog}\
        ${cmd} ${strain}
done

