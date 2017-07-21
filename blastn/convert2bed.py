import pandas as pd
import math

def main(allFilepath, strainFilepath):
    all_df=pd.read_csv(allFilepath)
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]

    for strain in strain_lst:
        strain_df=all_df[all_df["strain"]==strain]
        dct_lst=[]
        for _,row in strain_df.iterrows():
            dct={}
            dct["chrom"]=row["sseqid"]
            dct["chromStart"]=min(row["sstart"], row["send"])-1#-1 to 0 origin [start,end) format
            dct["chromEnd"]=max(row["sstart"], row["send"])

            dct["name"]="{0}:{1}".format(row["qseqid"].split(':')[0], row["family"])
            dct["score"]=int(-math.log(max(row["evalue"], 1e-10), 10)*100)
            dct["strand"]="."
            dct["thickStart"]=dct["chromStart"]
            dct["thickEnd"]=dct["chromEnd"]
            dct["itemRgb"]="255,0,0"
            dct_lst.append(dct)
        
        out_df=pd.DataFrame(dct_lst)
        out_df=out_df[["chrom","chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb"]]

        outFilepath="./bed/{0}.bed".format(strain)
        with open(outFilepath, 'w') as f:
            header="track name={0}_blastn.bed useScore=1".format(strain)
            f.write(header+'\n')
        out_df.to_csv(outFilepath, sep='\t', index=False, header=None, mode='a')

    
if __name__=="__main__":
    allFilepath="./out/all.csv"
    strainFilepath="/home/mitsuki/altorf/mbgd/data/strain.lst"
    main(allFilepath, strainFilepath)
