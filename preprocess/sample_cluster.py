import pandas as pd
import numpy

def get_ecoli_cluster(clusterFilepath, strain_lst=None):
    """
    get DataFrame of ecoli_cluster.tab
    when strain_lst is given, return relevant part
    """
    
    df=pd.read_csv(clusterFilepath,delimiter='\t', dtype="object")

    #format
    int_lst=["ClusterID", "HomClusterID","Size"]
    df[int_lst]=df[int_lst].astype(int)
    float_lst=["FuncCat(mbgd)", "FuncCat(cog)", "FuncCat(kegg)", "FuncCat(tigr)"]
    df[float_lst]=df[float_lst].astype(float)
    if strain_lst is None:
        strain_lst=list(df.columns[3:-7])
        assert len(strain_lst)==129

    #filter columns and add lineage column
    columns_lst=["family"]+strain_lst+["Description"]
    df=df[columns_lst]
    df["lineage"]=len(strain_lst)-df.isnull().sum(axis=1)
    
    return df[df["lineage"]>0]

def main(clusterFilepath, strainFilepath, sampledFilepath):
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    cluster_df=get_ecoli_cluster(clusterFilepath, strain_lst)
    
    dct_lst=[]
    for i in range(1,len(strain_lst)):
        filtered_df=cluster_df[cluster_df["lineage"]==i].sample(n=10, random_state=42)
        for key,row in filtered_df.iterrows():
            dct=dict(row)
            
            numQuery=0
            for strain in strain_lst:
                if isinstance(row[strain],str):
                    numQuery+=len(row[strain].split())
            dct["num_query"]=numQuery
            dct_lst.append(dct)

            
    out_df=pd.DataFrame(dct_lst)
    out_df=out_df[list(cluster_df.columns)+["num_query"]]
    out_df.to_csv(sampledFilepath,index=False)
    

if __name__=="__main__":
    clusterFilepath="/home/mitsuki/altorf/mbgd/data/ecoli_cluster.tab.mr"
    strainFilepath="/home/mitsuki/altorf/mbgd/data/strain.lst"
    sampledFilepath="/home/mitsuki/altorf/mbgd/data/sampled_cluster.csv"
    main(clusterFilepath, strainFilepath, sampledFilepath)
