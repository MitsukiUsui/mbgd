import pandas as pd


def get_ecoli_cluster(species_lst=None):
    """
    get DataFrame of ecoli_cluster.tab
    when species_lst is given, return relevant part
    """

    filepath="ecoli_cluster.tab"
    df=pd.read_csv(filepath,delimiter='\t', dtype="object")

    #format
    int_lst=["ClusterID", "HomClusterID","Size"]
    df[int_lst]=df[int_lst].astype(int)
    float_lst=["FuncCat(mbgd)", "FuncCat(cog)", "FuncCat(kegg)", "FuncCat(tigr)"]
    df[float_lst]=df[float_lst].astype(float)
    if species_lst is None:
        species_lst=list(df.columns[3:-7])
        assert len(species_lst)==129

    #filter columns and add lineage column
    columns_lst=["Family"]+species_lst+["Description"]
    df=df[columns_lst]
    df["lineage"]=len(species_lst)-df.isnull().sum(axis=1)
    
    return df[df["lineage"]>0]

def main():
    species_lst=[s.strip() for s in open("species.list", 'r').readlines()]
    cluster_df=get_ecoli_cluster(species_lst)
    
    dct_lst=[]
    for i in range(1,len(species_lst)):
        filtered_df=cluster_df[cluster_df["lineage"]==i].sample(n=10)
        for key,row in filtered_df.iterrows():
            dct=dict(row)
            
            numQuery=0
            for species in species_lst:
                if isinstance(row[species],str):
                    numQuery+=len(row[species].split())
            dct["num_query"]=numQuery
            dct_lst.append(dct)

            
    out_df=pd.DataFrame(dct_lst)
    out_df=out_df[list(cluster_df.columns)+["num_query"]]
    out_df.to_csv("sampled_cluster.csv",index=False)
    

if __name__=="__main__":
    main()