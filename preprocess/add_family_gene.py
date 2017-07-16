import pandas as pd
import re

def main(strainFilepath):
	
	#load cluster_df	
	clusterFilepath="ecoli_cluster.tab"
	cluster_df=pd.read_csv(clusterFilepath,delimiter='\t', dtype="object")
	int_lst=["ClusterID", "HomClusterID","Size"]
	cluster_df[int_lst]=cluster_df[int_lst].astype(int)
	float_lst=["FuncCat(mbgd)", "FuncCat(cog)", "FuncCat(kegg)", "FuncCat(tigr)"]
	cluster_df[float_lst]=cluster_df[float_lst].astype(float)

	strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
	for strain in strain_lst:
		print("PROCESSING: {}".format(strain))
		filepath="/data/mitsuki/data/mbgd/gene/{}.gene".format(strain)
		gene_df=pd.read_csv(filepath, delimiter='\t', header=None)
		
		gene_lst=list(set(gene_df[1]))
		
		#initialize dictionary with empty set
		gene2family_dct={}#map gene to set of family
		for gene in gene_lst:
			gene2family_dct[gene]=set()

		#update corresponding set
		pattern = r"([^()]+)(\([0-9]+\))?"
		r=re.compile(pattern)
		for _,row in cluster_df[["Family", strain]].iterrows():
			if isinstance(row[strain], str):
				for orfId in row[strain].split():
					gene=r.findall(orfId)[0][0].split(':')[1]#drop (num) and genome name
					gene2family_dct[gene].update([row["Family"]])#outer [] to add as string			  

		#create lookup_df from gene2family_dct
		dct_lst=[]
		for gene in gene_lst:
			dct={}
			dct[1]=gene
			if len(gene2family_dct[gene])==0:
				dct["Family"]="-1"
			else:
				dct["Family"]=','.join(list(gene2family_dct[gene]))
			dct_lst.append(dct)
		lookup_df=pd.DataFrame(dct_lst)

		gene_df=pd.merge(gene_df[gene_df.columns[:16]], lookup_df, on=1, how="left")
		gene_df.to_csv(filepath, index=False,header=None,sep='\t')


if __name__=="__main__":
	main("strain.lst")
	
