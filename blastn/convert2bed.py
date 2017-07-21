import pandas as pd
import math

def main(clusterFilepath, strainFilepath):
    cluster_df=pd.read_csv(clusterFilepath, delimiter='\t')
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]

    family_lst=list(cluster_df["Family"])
    
    for strain in strain_lst:
        print("PROCESSING {}".format(strain))
        dct_lst=[]
        for family in family_lst:
            filepath="/home/mitsuki/altorf/mbgd/blastn/result/{0}_{1}.tab".format(strain, family)
            if len(open(filepath, 'r').readlines())>0:
                df=pd.read_csv(filepath,delimiter='\t',header=None) 
                for _,row in df.iterrows():
                    dct={}
                    dct["chrom"]=row[1]
                    dct["chromStart"]=min(row[8], row[9])
                    dct["chromEnd"]=max(row[8], row[9])

                    dct["name"]="{0}:{1}".format(row[0].split(':')[0], family)
                    dct["score"]=int(-math.log(max(row[10], 1e-10), 10)*100)
                    dct["strand"]="."
                    dct["thickStart"]=dct["chromStart"]
                    dct["thickEnd"]=dct["chromEnd"]
                    dct["itemRgb"]="255,0,0"
                    dct_lst.append(dct)

        out_df=pd.DataFrame(dct_lst)
        out_df=out_df[["chrom","chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb"]]

        outFilepath="./out/{0}.bed".format(strain)
        with open(outFilepath, 'w') as f:
            header="track name={0}_blastn.bed useScore=1".format(strain)
            f.write(header+'\n')
        out_df.to_csv(outFilepath, sep='\t', index=False, header=None, mode='a')
    
if __name__=="__main__":
    clusterFilepath="/home/mitsuki/altorf/mbgd/data/sampled_cluster.tab"
    strainFilepath="/home/mitsuki/altorf/mbgd/data/strain.lst"
    main(clusterFilepath, strainFilepath)