from ete3 import Tree
import pandas as pd
import numpy as np

def main(clusterFilepath, strainFilepath, phbFilepath, outFilepath):
    cluster_df=pd.read_csv(clusterFilepath, dtype="object")
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    t=Tree(phbFilepath)
  
    print("Calc distance matrix")
    distance_mat=-np.ones((len(strain_lst), len(strain_lst)))
    for i, node1 in enumerate(strain_lst):
        for j, node2 in enumerate(strain_lst):
            if i!=j:
                distance_mat[i,j]=t.get_distance(node1, node2)

            
    dct_lst=[]
    for _, row in cluster_df.iterrows():
        if _%100==0:
            print(_)

        dct={}
        dct["family"]=row["family"]

        msk=row[strain_lst].isnull()
        for sidx in range(len(strain_lst)):
            if msk[sidx]:
                x = np.ma.array(distance_mat[sidx], mask=msk)
                qidx=x.argmin()
                assert distance_mat[sidx,qidx]>=0
                dct[strain_lst[sidx]]=strain_lst[qidx]
        dct_lst.append(dct)
    out_df=pd.DataFrame(dct_lst)
    out_df=out_df[["family"]+strain_lst]
    out_df.to_csv(outFilepath, index=False)


if __name__=="__main__":
    direc="../data/ecoli"
    clusterFilepath=direc+"/sampled_cluster.csv"
    strainFilepath=direc+"/strain.lst"
    phbFilepath=direc+"/tax562.phb"
    outFilepath="query.csv"
    main(clusterFilepath, strainFilepath, phbFilepath, outFilepath)