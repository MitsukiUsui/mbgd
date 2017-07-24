import pandas as pd
import math

def main():
    allFilepath="../blastn/out/all.csv"
    all_df=pd.read_csv(allFilepath)
    strainFilepath="/home/mitsuki/altorf/mbgd/data/strain.lst"
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    
    dct_lst=[]
    for strain in strain_lst: 
        geneFilepath="/data/mitsuki/data/mbgd/gene/{}.gene".format(strain)
        gene_df=pd.read_csv(geneFilepath, delimiter='\t', header=None)
        
        strain_df=all_df[all_df["sstrain"]==strain]#tmp_df for get chrname_lst
        chrname_lst=list(set(strain_df["sseqid"]))
        for chrname in chrname_lst:
            print("PROCESSING {}".format(chrname))
            filtered_df=strain_df[strain_df["sseqid"]==chrname]
        
            ###
            # 1. create ordered position of list together with CDS and blastn hits(regions).
            #    
            #    format: (pos, isStart, type, name)
            #    * sorted by first element(pos)
            #    * type: 0=CDS, 1=region
            ###
            pos_lst=[]
            for _, row in gene_df.iterrows():
                if row[7]=="CDS" and row[15]==chrname.split(':')[1]:
                    pos_lst.append((row[3], True,  0, row[1]))
                    pos_lst.append((row[4], False, 0, row[1]))
            for _, row in filtered_df.iterrows():
                pos_lst.append((min(row["sstart"],row["send"]), True,  1, row["region_id"]))
                pos_lst.append((max(row["sstart"],row["send"]),   False, 1, row["region_id"]))

            pos_lst=sorted(pos_lst, key=lambda x: (x[0], x[2]))
            
            
            ###
            # 2. find overlap simultaniously while scanning pos_lst in ascending order
            #    
            #    overlap_dctdct: key=(regionId, cdsname), value=dictionary of overlap information(ofirst, olast)
            #    cds_lst:    processing CDS
            #    region_lst: processing region
            ###
            overlap_dctdct={}
            cds_lst=[]
            region_lst=[]
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

            
            ###
            # 3. update dct_lst according to overlap_dctdct
            ###

            for k,v in overlap_dctdct.items():
                dct={}
                dct["region_id"]=k[0]
                dct["cds_name"]=k[1]
                dct["chr_name"]=chrname
                dct["ofirst"]=v["ofirst"]
                dct["olast"]=v["olast"]

                #calc sbjct first and last pos (gene, protein)
                row=gene_df[gene_df[1]==dct["cds_name"]].iloc[0,:]
                dct["cds_strand"]=row[5]
                if dct["cds_strand"]==1:
                    dct["sstart_gen"]=dct["ofirst"]-row[3]
                    dct["send_gen"]=dct["olast"]-row[3]+1
                else:
                    dct["sstart_gen"]=row[4]-dct["olast"]
                    dct["send_gen"]=row[4]-dct["ofirst"]+1
                dct["sstart_pro"]=math.ceil(dct["sstart_gen"]/3)
                dct["send_pro"]=math.floor(dct["send_gen"]/3)

                dct_lst.append(dct)
            print("\tfound {} overlaps".format(len(overlap_dctdct)))
                
    overlap_df=pd.DataFrame(dct_lst)
    overlap_df=overlap_df[["region_id", "chr_name", "ofirst", "olast", "cds_name","cds_strand", "sstart_gen", "send_gen", "sstart_pro", "send_pro"]]
    overlapFilepath="./out/overlap.csv"
    overlap_df.to_csv(overlapFilepath, index=False)
    print("DONE: {} overlaps in {}".format(overlap_df.shape[0], overlapFilepath))

    
if __name__=="__main__":
    main()