#!/home/mitsuki/.pyenv/versions/anaconda3-4.3.1/bin/python

import math
import sys
import numpy as np
import pandas as pd

def get_overlap_dctdct(gene_df, hit_df, chrname):
    """
    return overlap_dctdct

    key: (regionId, cdsname)
    value: dictionary of overlap information (ofirst, olast)
    """
    ###
    # 1. create ordered position of list together with CDS and blastn hits(regions).
    #    
    #    format: (pos, isStart, type, name)
    #    * sorted by first element(pos)
    #    * type: 0=CDS, 1=region
    ###
    filtered_df=hit_df[hit_df["sseqid"]==chrname]
    pos_lst=[]
    for _, row in gene_df.iterrows():
        if row[7]=="CDS" and row[15]==chrname.split(':')[1]:
            pos_lst.append((row[3], True,  0, row[1]))
            pos_lst.append((row[4], False, 0, row[1]))
    for _, row in filtered_df.iterrows():
        pos_lst.append((min(row["sstart"],row["send"]), True,  1, row["region_id"]))
        pos_lst.append((max(row["sstart"],row["send"]),   False, 1, row["region_id"]))
    pos_lst=sorted(pos_lst, key=lambda x: (x[0], x[2]))
    
    
    #find overlap simultaniously while scanning pos_lst in ascending order
    overlap_dctdct={}
    cds_lst=[]  #processing CDS
    region_lst=[]  #processing region
    for pos in pos_lst:
        #print(pos)
        if pos[2]==0:#if CDS
            if pos[1]:#if start
                cds_lst.append(pos[3])
                for region in region_lst:#add new overlap object with regions in region_lst
                    key=(region, pos[3])
                    overlap_dctdct[key]={}
                    overlap_dctdct[key]["ofirst"]=pos[0]
            else:#if last
                cds_lst.remove(pos[3])
                for region in region_lst:
                    key=(region, pos[3])
                    overlap_dctdct[key]["olast"]=pos[0]
        elif pos[2]==1:#if region
            if pos[1]:#if start
                region_lst.append(pos[3])
                for cds in cds_lst:#add new overlap object with CDSs in cds_lst
                    key=(pos[3], cds)
                    overlap_dctdct[key]={}
                    overlap_dctdct[key]["ofirst"]=pos[0]
            else:
                region_lst.remove(pos[3])
                for cds in cds_lst:
                    key=(pos[3], cds)
                    overlap_dctdct[key]["olast"]=pos[0]
    return overlap_dctdct

def add_sbjct_pos(overlap_df, gene_df):
    tmp_df=pd.DataFrame(index=gene_df.index)
    tmp_df["cds_name"]=gene_df[1]
    tmp_df["cds_first"]=gene_df[3]
    tmp_df["cds_last"]=gene_df[4]
    tmp_df["cds_strand"]=gene_df[5]
    tmp_df["sfamily"]=gene_df[16]
    merged_df=pd.merge(overlap_df, tmp_df, on="cds_name")
    assert merged_df.shape[0]==overlap_df.shape[0]

    dct_lst=[] 
    for _, row in merged_df.iterrows():
        dct={}
        if row["cds_strand"]==1:
            dct["sstart_dna"]=row["ofirst"]-row["cds_first"]
            dct["send_dna"]=row["olast"]-row["cds_first"]+1
        else:
            dct["sstart_dna"]=row["cds_last"]-row["olast"]
            dct["send_dna"]=row["cds_last"]-row["ofirst"]+1
        dct["sstart_pro"]=math.ceil(dct["sstart_dna"]/3)
        dct["send_pro"]=math.floor(dct["send_dna"]/3)
        dct_lst.append(dct)
    tmp_df=pd.DataFrame(dct_lst, index=merged_df.index)
    merged_df=pd.concat([merged_df, tmp_df], axis=1)
    
    merged_df=merged_df[["overlap_id","sstart_dna", "send_dna", "sstart_pro", "send_pro", 
                         "cds_strand", "sfamily"]]
    overlap_df=pd.merge(overlap_df, merged_df, on="overlap_id")
    return overlap_df

def add_query_pos(overlap_df, hit_df):
    ###
    # 1. calcurate the relative position of overlapping region to blastn hit, denoted by start_per, end_per each
    ###
    merged_df=pd.merge(overlap_df, hit_df, on="region_id")
    assert merged_df.shape[0]==overlap_df.shape[0]
    
    dct_lst=[]
    for _,row in merged_df.iterrows():
        dct={}
        length=abs(row["send"]-row["sstart"])
        if row["hit_strand"]==1:
            dct["start_per"]=(row["ofirst"]-row["sstart"])/length
            dct["end_per"]  =(row["olast"]-row["sstart"])/length
        else:
            dct["start_per"]=(row["sstart"]-row["olast"])/length
            dct["end_per"]  =(row["sstart"]-row["ofirst"])/length
        dct_lst.append(dct)
    tmp_df=pd.DataFrame(dct_lst, index=merged_df.index)
    merged_df=pd.concat([merged_df, tmp_df], axis=1)
    
    ###
    # 2. calculate qfirst_pro, qlast_pro
    ###
    merged_df["qstart_dna"]=np.ceil(merged_df["qstart"]+(merged_df["qend"]-merged_df["qstart"])*merged_df["start_per"]).astype(int)
    merged_df["qend_dna"]  =np.floor(merged_df["qstart"]+(merged_df["qend"]-merged_df["qstart"])*merged_df["end_per"]).astype(int)
    merged_df["qstart_pro"]=np.ceil(merged_df["qstart_dna"]/3).astype(int)
    merged_df["qend_pro"]  =np.floor(merged_df["qend_dna"]/3).astype(int)

    merged_df["qorf_id"]=merged_df["qseqid"]
    merged_df=merged_df[["overlap_id","start_per", "end_per", "qstart_dna", "qend_dna", "qstart_pro", "qend_pro", 
                         "qorf_id", "qfamily", "hit_strand"]]
    overlap_df=pd.merge(overlap_df, merged_df, on="overlap_id")
    return overlap_df

def main(strain, hitFilepath, geneFilepath, overlapFilepath):
    hit_df=pd.read_csv(hitFilepath)
    gene_df=pd.read_csv(geneFilepath, delimiter='\t', header=None)

    chrname_lst=list(set(hit_df["sseqid"]))
    dct_lst=[]
    for chrname in chrname_lst:
        print("PROCESSING {}".format(chrname))
        overlap_dctdct=get_overlap_dctdct(gene_df, hit_df, chrname)
        
        #update dct_lst according to overlap_dctdct
        for k,v in overlap_dctdct.items():
            dct={}
            dct["region_id"]=k[0]
            dct["cds_name"]=k[1]
            dct["chr_name"]=chrname
            dct["ofirst"]=v["ofirst"]
            dct["olast"]=v["olast"]
            dct["olength"]=v["olast"]-v["ofirst"]+1
            dct_lst.append(dct)
        print("\tfound {} overlaps".format(len(overlap_dctdct)))
        
    overlap_df=pd.DataFrame(dct_lst)
    overlap_df["overlap_id"]=overlap_df.index
    overlap_df=add_sbjct_pos(overlap_df, gene_df) #"sstart_dna","send_dna", "sstart_pro", "send_pro" + alpha
    overlap_df=add_query_pos(overlap_df, hit_df) #"qstart_dna","qend_dna", "qstart_pro", "qend_pro", "start_per", "end_per" + alpha
    overlap_df["sorf_id"]=strain + ':' + overlap_df["cds_name"]
   
    columns_lst=["overlap_id","ofirst", "olast", "olength", "start_per", "end_per", 
                 "qfamily","sfamily","qorf_id", "sorf_id", "hit_strand", "cds_strand",
                 "qstart_dna","qend_dna", "qstart_pro", "qend_pro", 
                 "sstart_dna","send_dna", "sstart_pro", "send_pro",
                 "region_id", "chr_name"] 
    overlap_df=overlap_df[columns_lst]
    overlap_df=overlap_df.sort_values(by="overlap_id")
    overlap_df.to_csv(overlapFilepath, index=False)
    print("DONE: {} overlaps in {}".format(overlap_df.shape[0], overlapFilepath))


if __name__=="__main__":
    strain=sys.argv[1]
    hitFilepath="/home/mitsuki/altorf/mbgd/blastn/out/{}.csv".format(strain)
    geneFilepath="/data/mitsuki/data/mbgd/gene/{}.gene".format(strain)
    overlapFilepath="/home/mitsuki/altorf/mbgd/analyze/out/{}_ovr.csv".format(strain)
    main(strain,hitFilepath, geneFilepath, overlapFilepath)
