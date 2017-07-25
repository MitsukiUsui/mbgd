import pandas as pd
import numpy as np
import subprocess

def check(clusterFilepath, strainFilepath):
	cluster_df=pd.read_csv(clusterFilepath, delimiter='\t', dtype="object")
	strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]

	familyCount=0
	jobCount=0
	resultCount=0

	for _,row in cluster_df.iterrows():
		msk=list(row[strain_lst].isnull())# True if the strain does not have the family
		size=len(strain_lst)-np.sum(msk)
		if 0<size and size<len(strain_lst):
			familyCount+=1
			for i, strain in enumerate(strain_lst):
				if msk[i]:
					jobCount+=1
					filepath="/home/mitsuki/altorf/mbgd/blastn/result/{}_{}.tab".format(strain, row["family"])
					try:
						open(filepath, 'r')
						resultCount+=1
					except FileNotFoundError:
						pass
	
	print("TOTAL")	
	print("\tfamily: {}".format(familyCount))
	print("\tjobs:   {}/{}".format(resultCount, jobCount))

def main(clusterFilepath, strainFilepath):
	cluster_df=pd.read_csv(clusterFilepath, delimiter='\t', dtype="object")
	strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]

	for _,row in cluster_df.iterrows():
		msk=list(row[strain_lst].isnull())# True if the strain does not have the family
		size=len(strain_lst)-np.sum(msk)
		if 0<size and size<len(strain_lst):
			for i, strain in enumerate(strain_lst):
				if msk[i]:
					cmd="/home/mitsuki/altorf/mbgd/blastn/uge_blastn_args.sh {} {}".format(strain,row["family"])
					subprocess.check_call(cmd.strip().split(' '))

if __name__=="__main__":
	#clusterFilepath="/home/mitsuki/altorf/mbgd/data/sampled_cluster.tab"
	clusterFilepath="/home/mitsuki/altorf/mbgd/data/ecoli_cluster.tab.mr"
	strainFilepath="/home/mitsuki/altorf/mbgd/data/strain.lst"
	main(clusterFilepath, strainFilepath)
	#check(clusterFilepath, strainFilepath)
			
