import pandas as pd

def main(filepath):
    cluster_df=pd.read_csv(filepath, delimiter='\t', dtype="object")
    cluster_df["family"]=["family{}".format(i) for i in cluster_df["ClusterID"]]
    assert len(set(cluster_df["family"]))==cluster_df.shape[0]

    cluster_df["alias"]=cluster_df["Gene"]
    cluster_df["alias"]=cluster_df["alias"].fillna(cluster_df["family"])
    cluster_df.to_csv(filepath, index=False, sep='\t')
    print("OUTPUT to {}".format(filepath))

if __name__=="__main__":
    filepath="~/altorf/mbgd/data/streptomyces/streptomyces_cluster.tab"
    main(filepath)
