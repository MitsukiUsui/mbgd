strain=${1}
dir=/home/mitsuki/altorf/mbgd/blastn
dbName=${dir}/db/${strain}
queryFilepath=${dir}/query/${strain}.query
outFilepath=${dir}/result/${strain}.tab

blastn -db ${dbName}\
	   -query ${queryFilepath}\
	   -out ${outFilepath}\
	   -word_size 6\
	   -evalue 1e-3\
	   -outfmt 6 

