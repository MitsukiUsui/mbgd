import pandas as pd
import re

if __name__=="__main__":
	filepath="ecoli_cluster.tab"
	
	cluster_df=pd.read_csv(filepath,delimiter='\t', dtype="object")
	int_lst=["ClusterID", "HomClusterID","Size"]
	cluster_df[int_lst]=cluster_df[int_lst].astype(int)
	float_lst=["FuncCat(mbgd)", "FuncCat(cog)", "FuncCat(kegg)", "FuncCat(tigr)"]
	cluster_df[float_lst]=cluster_df[float_lst].astype(float)


	lst=[]
	for key, row in cluster_df.iterrows():
		if isinstance(row["Gene"],str):
			lst.append(row["Gene"])
		else:
			lst.append("Family{}".format(row["ClusterID"]))
	cluster_df["Family"]=lst

	#filepath="ecoli_cluster.tab.new"
	cluster_df.to_csv(filepath, sep='\t', index=False)
