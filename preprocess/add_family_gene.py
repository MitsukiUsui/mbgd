#!/usr/bin/env python3

import pandas as pd
import re
import sys
import os

def main(target, clusterFilepath,strainFilepath):
    #load cluster_df    
    cluster_df=pd.read_csv(clusterFilepath, delimiter='\t', dtype="object")
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    inDirec="/data/mitsuki/data/mbgd/gene"
    outDirec="/data/mitsuki/data/mbgd/{}/gene".format(target)
    os.makedirs(outDirec, exist_ok=True)
    
    for strain in strain_lst:
        print("START: process {}".format(strain))
        
        inFilepath="{}/{}.gene".format(inDirec, strain)
        gene_df=pd.read_csv(inFilepath, delimiter='\t', header=None)
        gene_lst=list(set(gene_df[1]))
        
        #initialize dictionary with empty set
        gene2family_dct={}#map gene to set of family
        for gene in gene_lst:
            gene2family_dct[gene]=set()

        #update gene2family_dct, iterating each row
        pattern = r"([^()]+)(\([0-9]+\))?"
        r=re.compile(pattern)
        for _,row in cluster_df[["family", strain]].iterrows():
            if isinstance(row[strain], str):
                for orfId in row[strain].split():
                    gene=r.findall(orfId)[0][0].split(':')[1]#drop (num) and genome name
                    gene2family_dct[gene].update([row["family"]])#outer [] to add as string           

        #create lookup_df from gene2family_dct
        dct_lst=[]
        for gene in gene_lst:
            dct={}
            dct[1]=gene
            if len(gene2family_dct[gene])>0:
                dct["family"]=','.join(list(gene2family_dct[gene]))
            dct_lst.append(dct)
        lookup_df=pd.DataFrame(dct_lst)
        
        # join lookup_df to gene_df
        gene_df=pd.merge(gene_df, lookup_df, on=1, how="left")
        print("\tadded family to {}/{} genes".format((~gene_df["family"].isnull()).sum(),  gene_df.shape[0]))
        
        # output
        outFilepath="{}/{}.gene".format(outDirec, strain)
        gene_df.to_csv(outFilepath, index=False,header=None,sep='\t')
        print("\tDONE: output to {}".format(outFilepath))

if __name__=="__main__":
    target=sys.argv[1]
    direc="/home/mitsuki/altorf/mbgd/data/{}".format(target)
    clusterFilepath="{}/cluster.tsv".format(direc)
    strainFilepath="{}/strain.lst".format(direc)
    main(target, clusterFilepath, strainFilepath)
