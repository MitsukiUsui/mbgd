dir=/home/mitsuki/altorf/mbgd/blastn

sstrain=${1}
family=${2}
qstrain=${3}

dbName=${dir}/db/${sstrain}
queryFilepath=/data/mitsuki/data/mbgd/family/geneseq/${family}/${family}_${qstrain}.geneseq
outFilepath=${dir}/result/${sstrain}_${family}.tab

blastn -db ${dbName}\
	   -query ${queryFilepath}\
	   -out ${outFilepath}\
	   -word_size 6\
	   -evalue 1e-3\
	   -outfmt 6 

