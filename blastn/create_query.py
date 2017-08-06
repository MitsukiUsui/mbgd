import pandas as pd
from Bio import SeqIO

def main(lookupFilepath, strainFilepath):
    lookup_df=pd.read_csv(lookupFilepath, index_col=0)
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    
    for sstrain in strain_lst:
        seqRec_lst=[]
        for family, orfId in lookup_df[sstrain].dropna().iteritems():
            qstrain=orfId.split(':')[0]
            seqFilepath="/data/mitsuki/data/mbgd/family/geneseq/{0}/{0}_{1}.geneseq".format(family, qstrain)
            for seqRec in SeqIO.parse(seqFilepath, "fasta"):
                seqRec_lst.append(seqRec)
        
        queryFilepath="./query/{}.query".format(sstrain)
        with open(queryFilepath, 'w') as f:
            for seqRec in seqRec_lst:
                SeqIO.write(seqRec, f, "fasta")
        print("OUTPUT to {}".format(queryFilepath))


if __name__=="__main__":
    lookupFilepath="./out/query_lookup.csv"
    strainFilepath="../data/ecoli/strain.lst"
    main(lookupFilepath, strainFilepath)
