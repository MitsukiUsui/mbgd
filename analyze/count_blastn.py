import pandas as pd

def main():
    strainFilepath="../preprocess/strain.lst"
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    cluster_df=pd.read_csv("../preprocess/sampled_cluster.csv")
    family_lst=list(cluster_df["Family"])
    
    dct_lst=[]
    for family in family_lst:
        dct={}
        dct["Family"]=family
        for strain in strain_lst:
            filepath="../blastn/result/{0}_{1}.tab".format(strain,family)
            dct[strain]=len(open(filepath, 'r').readlines())
        dct_lst.append(dct)
        
    count_df=pd.DataFrame(dct_lst)
    count_df=count_df[["Family"]+strain_lst]
    count_df["lineage"]=cluster_df["lineage"]
    count_df["num_query"]=cluster_df["num_query"]
    count_df.to_csv("count_blastn.csv",index=False)


if __name__=="__main__":
    main()