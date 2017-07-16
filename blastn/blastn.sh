dir=/home/mitsuki/altorf/mbgd/blastn

strain=${1}
family=${2}

dbName=${dir}/db/${strain}
queryFilepath=${dir}/query/${family}.fna
outFilepath=${dir}/result/${strain}_${family}.tab

blastn -db ${dbName}\
	   -query ${queryFilepath}\
	   -out ${outFilepath}\
	   -word_size 6\
	   -evalue 1e-3\
	   -outfmt 6 

#outFilepath=${dir}/out/${strain}_${family}.blasttxt
#blastn -db ${dbName}\
#	   -query ${queryFilepath}\
#	   -out ${outFilepath}\
#	   -word_size 6\
#	   -evalue 1e-3 

