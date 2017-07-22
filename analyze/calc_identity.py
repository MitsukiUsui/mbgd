import pandas as pd
import numpy as np
import sys
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

sys.path.append("/home/mitsuki/altorf/mbgd/blastn")
from create_query import extract_sequence

def main():
    allFilepath="../blastn/out/all.csv"
    all_df=pd.read_csv(allFilepath)
    overlapFilepath="./out/overlap.csv"
    overlap_df=pd.read_csv(overlapFilepath)
    merged_df=pd.merge(overlap_df, all_df, how="left", on="region_id")
    
    merged_df["strand"]='+'
    merged_df.loc[(merged_df["sstart"]>merged_df["send"]), "strand"]='-'
    
    #calculate first_per, last_per    
    tpl_lst=[]
    for _,row in merged_df.iterrows():
        length=abs(row["send"]-row["sstart"])
        if row["strand"]=='+':
            first_per=(row["ofirst"]-row["sstart"])/length
            last_per=(row["olast"]-row["sstart"])/length
        else:
            first_per=(row["sstart"]-row["olast"])/length
            last_per=(row["sstart"]-row["ofirst"])/length
        tpl_lst.append((first_per,last_per))
    merged_df["first_per"]=0
    merged_df["last_per"]=0
    merged_df[["first_per","last_per"]]=tpl_lst
    
    #calculate qfirst_pro, qlast_pro
    merged_df["qfirst_pro"]=np.ceil((merged_df["qstart"]+(merged_df["qend"]-merged_df["qstart"])*merged_df["first_per"])/3).astype(int)
    merged_df["qlast_pro"]=np.floor((merged_df["qstart"]+(merged_df["qend"]-merged_df["qstart"])*merged_df["last_per"])/3).astype(int)
    
    #global alignment for each overlap
    matrix = matlist.blosum62
    gap_open=-10
    gap_extend=-0.5
    
    lst=[]
    for key,row in merged_df.iterrows():
        if key%100==0:
            print(key)
            
        if row["olast"]-row["ofirst"]>10:
            seq=extract_sequence(row["qseqid"],dtype="proteinseq")
            qseq=seq[max(0,row["qfirst_pro"]):min(len(seq)-1,row["qlast_pro"])]#-1 to avoid inclusion of '*'
            seq=extract_sequence(row["chr_name"].split(':')[0]+':'+row["cds_name"],dtype="proteinseq")
            sseq=seq[max(0,row["sfirst_pro"]):min(len(seq)-1,row["slast_pro"])]

            print(key)
            print(row["qfirst_pro"],row["qlast_pro"])
            print(qseq)
            print(row["sfirst_pro"],row["slast_pro"])
            print(sseq)
            print()

            alns=pairwise2.align.globalds(qseq, sseq, matrix, gap_open, gap_extend)
            score=alns[0][2]
            lst.append(score)
        else:
            lst.append(0)
    merged_df["score"]=lst
    merged_df.to_csv("./out/merged.csv",index=False)

if __name__=="__main__":
    main()
    