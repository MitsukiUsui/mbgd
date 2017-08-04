import pandas as pd
import re

def main(clusterFilepath,strainFilepath):
    #load cluster_df    
    cluster_df=pd.read_csv(clusterFilepath, delimiter='\t', dtype="object")

    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    for strain in strain_lst:
        print("PROCESSING {}".format(strain))
        filepath="/data/mitsuki/data/mbgd/gene/{}.gene".format(strain)
        gene_df=pd.read_csv(filepath, delimiter='\t', header=None)
        
        gene_lst=list(set(gene_df[1]))
        
        #initialize dictionary with empty set
        gene2family_dct={}#map gene to set of family
        for gene in gene_lst:
            gene2family_dct[gene]=set()

        #update gene2family_dct, iterating each row
        pattern = r"([^()]+)(\([0-9]+\))?"
        r=re.compile(pattern)
        for _,row in cluster_df[["family", strain]].iterrows():
            if isinstance(row[strain], str):
                for orfId in row[strain].split():
                    gene=r.findall(orfId)[0][0].split(':')[1]#drop (num) and genome name
                    gene2family_dct[gene].update([row["family"]])#outer [] to add as string           

        #create lookup_df from gene2family_dct
        dct_lst=[]
        count=0
        
        for gene in gene_lst:
            dct={}
            dct[1]=gene
            if len(gene2family_dct[gene])>0:
                dct["family"]=','.join(list(gene2family_dct[gene]))
                count+=1
            dct_lst.append(dct)
        lookup_df=pd.DataFrame(dct_lst)
        gene_df=pd.merge(gene_df[gene_df.columns[:16]], lookup_df, on=1, how="left")
        gene_df.to_csv(filepath, index=False,header=None,sep='\t')
        print("\tadded family to {}/{} genes".format(count, len(gene_lst)))
              

if __name__=="__main__":
    direc="../data/streptomyces"
    clusterFilepath=direc+"/streptomyces_cluster.tab"
    strainFilepath=direc+"/strain.lst"
    main(clusterFilepath, strainFilepath)
