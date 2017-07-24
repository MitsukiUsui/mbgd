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

def output_all(cluster_df, strain_lst, allFilepath):
	columns_lst=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
							 "qstart", "qend", "sstart","send", "evalue", "bitscore"]
	
	dct_lst=[]
	for _,row in cluster_df.iterrows():
		msk=list(row[strain_lst].isnull())# True if the strain does not have the family
		for i, strain in enumerate(strain_lst):
			if msk[i]:
				filepath="./result/{0}_{1}.tab".format(strain,row["family"])
				if len(open(filepath, 'r').readlines())>0:
					df=pd.read_csv(filepath,delimiter='\t', header=None)
					for _,subrow in df.iterrows():
						dct={}
						dct["qfamily"]=row["family"]
						dct["qstrain"]=subrow[0].split(':')[0]
						dct["sstrain"]=subrow[1].split(':')[0]

						for i,column in enumerate(columns_lst):
							dct[column]=subrow[i]
						dct_lst.append(dct)
	
	all_df=pd.DataFrame(dct_lst)
	all_df=all_df[["qfamily", "qstrain", "sstrain"]+columns_lst]
	all_df.index.name="region_id"
	all_df.to_csv(allFilepath)
	

if __name__=="__main__":
	clusterFilepath="/home/mitsuki/altorf/mbgd/data/sampled_cluster.tab"
	strainFilepath="/home/mitsuki/altorf/mbgd/data/strain.lst"
   
	strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
	cluster_df=pd.read_csv(clusterFilepath, delimiter='\t')

	countFilepath="./out/count.csv"
	output_count(cluster_df, strain_lst, countFilepath)
	
	allFilepath="./out/all.csv"
	output_all(cluster_df, strain_lst, allFilepath)
