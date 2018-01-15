target=${1}
dataDirec=/home/mitsuki/altorf/mbgd/data/${target}
outFilepath=${dataDirec}/cluster.phb

#--------------------------------------------------------------------------------
# force mode
#--------------------------------------------------------------------------------
FORCE_MODE=false
forceFilepath=${outFilepath}
if [ "$FORCE_MODE" = false ] && [ -e ${forceFilepath} ]; then
    echo "PASS: target file already exists"
    exit
fi

phyloDirec=/home/mitsuki/software/phylophlan
mkdir -p ${phyloDirec}/input/${target}
strainFilepath=${dataDirec}/strain.lst
while read strain
do
    from=/data/mitsuki/data/mbgd/proteinseq/${strain}.proteinseq
    to=${phyloDirec}/input/${target}/${strain}.faa
    ln -s ${from} ${to}
done < ${strainFilepath}

cd ${phyloDirec}
./phylophlan.py -u ${target} --nproc 4

nwkFilepath=${phyloDirec}/output/${target}/${target}.tree.nwk
ln -s ${nwkFilepath} ${outFilepath}
