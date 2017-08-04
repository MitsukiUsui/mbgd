import pandas as pd

def output_count(cluster_df, strain_lst, countFilepath):
    dct_lst=[]
    for _,row in cluster_df.iterrows():
        dct={}
        dct["family"]=row["family"]
        dct["lineage"]=row["lineage"]
        dct["num_query"]=row["num_query"]
        
        msk=list(row[strain_lst].isnull())# True if the strain does not have the family
        for i, strain in enumerate(strain_lst):
            if msk[i]:
                filepath="./result/{0}_{1}.tab".format(strain,row["family"])
                dct[strain]=len(open(filepath, 'r').readlines())
            else:
                dct[strain]=-1 
        dct_lst.append(dct)
    count_df=pd.DataFrame(dct_lst)
    count_df=count_df[["family"]+strain_lst+["lineage", "num_query"]]
    count_df.to_csv(countFilepath,index=False)
    print("DONE counting to {}".format(countFilepath))

def output_strain(strain_lst, family_lst):
    columns_lst=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                             "qstart", "qend", "sstart","send", "evalue", "bitscore"]

    for strain in strain_lst:
        strainFilepath="./out/{}.csv".format(strain)
        dct_lst=[]

        for family in family_lst:
            filepath="./result/{0}_{1}.tab".format(strain,family)
            try:
                if len(open(filepath, 'r').readlines())>0:
                    df=pd.read_csv(filepath,delimiter='\t', header=None)
                    for _,row in df.iterrows():
                        dct={}
                        dct["qfamily"]=family
                        dct["qstrain"]=row[0].split(':')[0]
                        dct["sstrain"]=row[1].split(':')[0]
                        if row[8]<row[9]:
                            dct["hit_strand"]=1
                        else:
                            dct["hit_strand"]=-1
                        
                        for i,column in enumerate(columns_lst):
                            dct[column]=row[i]
                        dct_lst.append(dct)
            except FileNotFoundError:
                pass
        
        strain_df=pd.DataFrame(dct_lst)
        strain_df=strain_df[["qfamily", "qstrain", "sstrain"]+columns_lst+["hit_strand"]]
        strain_df.index.name="region_id"
        strain_df.to_csv(strainFilepath)
        print("DONE integration to {}".format(strainFilepath))
    

if __name__=="__main__":
    direc="../data/streptomyces"
    clusterFilepath=direc+"/sampled_cluster.csv"
    strainFilepath=direc+"/strain.lst"
   
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    cluster_df=pd.read_csv(clusterFilepath, dtype="object")

    countFilepath="./out/test_count.csv"
    output_count(cluster_df, strain_lst, countFilepath)

    family_lst=list(cluster_df["family"])
    output_strain(strain_lst,family_lst)
