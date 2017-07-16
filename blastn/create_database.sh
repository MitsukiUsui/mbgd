IFS=$'\n'

for strain in `cat ../preprocess/strain.list`
do
	inFilepath=/data/mitsuki/data/mbgd/dnaseq/${strain}.dnaseq
	dbName=~/altorf/mbgd/blastn/db/${strain}

	makeblastdb -in ${inFilepath} -dbtype nucl -out ${dbName}
	echo "CREATED ${dbName}"
done
