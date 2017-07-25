import pandas as pd
import re
from collections import Counter

if __name__=="__main__":
	filepath="~/altorf/mbgd/data/ecoli_cluster.tab"
	
	cluster_df=pd.read_csv(filepath,delimiter='\t', dtype="object")
	int_lst=["ClusterID", "HomClusterID","Size"]
	cluster_df[int_lst]=cluster_df[int_lst].astype(int)
	float_lst=["FuncCat(mbgd)", "FuncCat(cog)", "FuncCat(kegg)", "FuncCat(tigr)"]
	cluster_df[float_lst]=cluster_df[float_lst].astype(float)

	#count #occurence of each cluster
	family_lst=list(cluster_df["Gene"].dropna())
	counter_dct=Counter(family_lst)
	for k,v in counter_dct.items():
		if v==1:
			counter_dct[k]=-1
	
	lst=[]
	for _, row in cluster_df.iterrows():
		if isinstance(row["Gene"],str):
			name=row["Gene"]
			if counter_dct[name]==-1:
				lst.append(name)
			else:
				lst.append("{}#{}".format(name,counter_dct[name]))
				counter_dct[name]-=1
		else:
			lst.append("family{}".format(row["ClusterID"]))
	lst=[s.replace('/','_') for s in lst]
	cluster_df["family"]=lst

	outFilepath="~/altorf/mbgd/data/ecoli_cluster.tab.mr"
	cluster_df.to_csv(outFilepath, sep='\t', index=False)
