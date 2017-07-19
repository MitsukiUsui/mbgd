import pandas as pd
import subprocess

def main(clusterFilepath, strainFilepath):
    cluster_df=pd.read_csv(clusterFilepath, delimiter='\t')
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]

    for _,row in cluster_df.iterrows():
        msk=list(row[strain_lst].isnull())# True if the strain does not have the family
        for i, strain in enumerate(strain_lst):
            if msk[i]:
                cmd="/home/mitsuki/altorf/mbgd/blastn/sge_blastn_args.sh {} {}".format(strain,row["family"])
                print(cmd)
                subprocess.check_call(cmd.strip().split(' '))
            

if __name__=="__main__":
    clusterFilepath="/home/mitsuki/altorf/mbgd/data/sampled_cluster.tab"
    strainFilepath="/home/mitsuki/altorf/mbgd/data/strain.lst"
    main(clusterFilepath, strainFilepath)
            
