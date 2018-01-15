#!/usr/bin/env python3

"""
add 3 columns below to cluster.tab and output to cluster.tsv
    family...family id 
    alias...gene name or family id
    lineage...#strain which has the family
"""

import pandas as pd
import sys

def main(strainFilepath, inFilepath, outFilepath):
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    
    print("START: load {}".format(inFilepath))
    cluster_df=pd.read_csv(inFilepath, delimiter='\t', dtype="object")
    
    cluster_df["family"]=["family{}".format(i) for i in cluster_df["ClusterID"]]
    assert len(set(cluster_df["family"]))==cluster_df.shape[0]
    cluster_df["alias"]=cluster_df["Gene"]
    cluster_df["alias"]=cluster_df["alias"].fillna(cluster_df["family"])
    cluster_df["lineage"]=len(strain_lst)-cluster_df[strain_lst].isnull().sum(axis=1)
    
    print("DONE: add 3 new columns (family, alias, lineage)")
    
    cluster_df.to_csv(outFilepath, index=False, sep='\t')
    print("DONE: output to {}".format(outFilepath))
    
if __name__=="__main__":
    target=sys.argv[1]
    direc="../data/{}".format(target)
    strainFilepath="{}/strain.lst".format(direc)
    inFilepath="{}/cluster.tab".format(direc)
    outFilepath="{}/cluster.tsv".format(direc)
    main(strainFilepath, inFilepath, outFilepath)
    