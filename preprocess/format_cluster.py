import pandas as pd
import numpy

def read_ecoli_cluster(clusterFilepath):
	cluster_df=pd.read_csv(clusterFilepath,delimiter='\t', dtype="object")
	int_lst=["ClusterID", "HomClusterID","Size"]
	cluster_df[int_lst]=cluster_df[int_lst].astype(int)
	float_lst=["FuncCat(mbgd)", "FuncCat(cog)", "FuncCat(kegg)", "FuncCat(tigr)"]
	cluster_df[float_lst]=cluster_df[float_lst].astype(float)
	return cluster_df 

def format_cluster(cluster_df, strain_lst=None):
	"""
	adds lineage, num_query column
	deletes unnecessary columns and rows
	"""

	if strain_lst is None:
		strain_lst=list(cluster_df.columns[3:-7])
		assert len(strain_lst)==129

	cluster_df["lineage"]=len(strain_lst)-cluster_df[strain_lst].isnull().sum(axis=1)
	cluster_df=cluster_df[cluster_df["lineage"]>0].copy()#delete unrelated rows

	#add num_query column
	lst=[]
	for _,row in cluster_df.iterrows():
		numQuery=0
		for strain in strain_lst:
			if isinstance(row[strain],str):
				numQuery+=len(row[strain].split())
		lst.append(numQuery)
	cluster_df["num_query"]=lst
	
	columns_lst=["family"]+strain_lst+["Description", "lineage", "num_query"]
	cluster_df=cluster_df[columns_lst]
	return cluster_df


def sample_cluster(cluster_df, numStrain, n=10):
	df=pd.DataFrame(columns=cluster_df.columns, dtype=cluster_df.dtypes)
	for i in range(1,numStrain):
		filtered_df=cluster_df[cluster_df["lineage"]==i].sample(n=n, random_state=42)
		df=pd.concat([df, filtered_df])
	return df

def main():
	#configuration
	clusterFilepath="/home/mitsuki/altorf/mbgd/data/ecoli_cluster.tab.mr"
	strainFilepath="/home/mitsuki/altorf/mbgd/data/strain.lst"
	outFilepath="/home/mitsuki/altorf/mbgd/data/test_cluster.csv"
	sample=True
	
	cluster_df=read_ecoli_cluster(clusterFilepath)
	strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
	cluster_df=format_cluster(cluster_df, strain_lst)
	
	if sample:
		cluster_df=sample_cluster(cluster_df, len(strain_lst))	
	
	cluster_df["lineage"]=cluster_df["lineage"].astype(int)
	cluster_df["num_query"]=cluster_df["num_query"].astype(int)
	cluster_df.to_csv(outFilepath, index=False)
	print("DONE output to {}".format(outFilepath))

if __name__=="__main__":
	main()
