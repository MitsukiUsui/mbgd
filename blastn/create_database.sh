#!/bin/bash

IFS=$'\n'

target=${1}
strainFilepath=../data/${target}/strain.lst
outDirec=./db/${target}
mkdir -p ${outDirec}

for strain in `cat ${strainFilepath}`
do
	inFilepath=/data/mitsuki/data/mbgd/dnaseq/${strain}.dnaseq
	dbName=${outDirec}/${strain}
	makeblastdb -in ${inFilepath} -dbtype nucl -out ${dbName}
	echo "CREATED ${dbName}"
done
