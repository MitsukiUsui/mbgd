import pandas as pd

def main():
    allFilepath="../blastn/out/all.csv"
    all_df=pd.read_csv(allFilepath)
    strainFilepath="/home/mitsuki/altorf/mbgd/data/strain.lst"
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    
    dct_lst=[]
    for strain in strain_lst: 
        geneFilepath="/data/mitsuki/data/mbgd/gene/{}.gene".format(strain)
        gene_df=pd.read_csv(geneFilepath, delimiter='\t', header=None)
        
        strain_df=all_df[all_df["strain"]==strain]#tmp_df for get chrname_lst
        chrname_lst=list(set(strain_df["sseqid"]))
        for chrname in chrname_lst:
            print("PROCESSING {}".format(chrname))
            
            filtered_df=strain_df[strain_df["sseqid"]==chrname]
        
            pos_lst=[]#position list (pos, isStart, type, name)
            for _, row in gene_df.iterrows():
                if row[7]=="CDS" and row[15]==chrname.split(':')[1]:
                    pos_lst.append((row[3], True,  0, row[1]))
                    pos_lst.append((row[4], False, 0, row[1]))
            for _, row in filtered_df.iterrows():
                pos_lst.append((min(row["sstart"],row["send"]), True,  1, row["region_id"]))
                pos_lst.append((max(row["sstart"],row["send"]),   False, 1, row["region_id"]))
            
            
            pos_lst=sorted(pos_lst, key=lambda x: (x[0], x[2]))
            
            
            #update overlap_dctdct by processing pos_lst in ascending order
            overlap_dctdct={}#key:(regionId, cdsname) value:{first, last}
            cds_lst=[]#remember processing CDS
            region_lst=[]#                 region
            for pos in pos_lst:#t for type
                #print(pos)
                if pos[2]==0:#if CDS
                    if pos[1]:#if start
                        cds_lst.append(pos[3])
                        for region in region_lst:
                            key=(region, pos[3])
                            overlap_dctdct[key]={}
                            overlap_dctdct[key]["first"]=pos[0]
                    else:#if last
                        cds_lst.remove(pos[3])
                        for region in region_lst:
                            key=(region, pos[3])
                            overlap_dctdct[key]["last"]=pos[0]
                elif pos[2]==1:#if region
                    if pos[1]:#if start
                        region_lst.append(pos[3])
                        for cds in cds_lst:
                            key=(pos[3], cds)
                            overlap_dctdct[key]={}
                            overlap_dctdct[key]["first"]=pos[0]
                    else:
                        region_lst.remove(pos[3])
                        for cds in cds_lst:
                            key=(pos[3], cds)
                            overlap_dctdct[key]["last"]=pos[0]

            #update dct_lst according to overlap_dctdct
            for k,v in overlap_dctdct.items():
                dct={}
                dct["region_id"]=k[0]
                dct["chr_name"]=chrname
                dct["cds_name"]=k[1]
                dct["first"]=v["first"]
                dct["last"]=v["last"]
                dct_lst.append(dct)
            print("\tfound {} overlaps".format(len(overlap_dctdct)))
                
    overlap_df=pd.DataFrame(dct_lst)
    overlap_df=overlap_df[["region_id", "chr_name", "cds_name", "first", "last"]]
    overlapFilepath="overlap.csv"
    overlap_df.to_csv(overlapFilepath, index=False)
    print("DONE: {} overlaps in {}".format(overlap_df.shape[0], overlapFilepath))

    
if __name__=="__main__":
    main()