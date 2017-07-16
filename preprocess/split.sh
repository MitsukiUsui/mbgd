dir=/data/mitsuki/data/mbgd
inFilepath=${dir}/mbgd_2016-01_gene

awk -F "\t" -v d="${dir}" 'NR>1{fn=sprintf("%s/gene/%s.gene",d,$1); print $0 >> fn}' ${inFilepath}
