import pandas as pd

def main(lookupFilepath, strainFilepath):
    lookup_df=pd.read_csv(lookupFilepath)
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()] 
    for sstrain in strain_lst:
        hitFilepath="./result/{}.tab".format(sstrain)
        columns_lst=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", 
                     "sstart","send", "evalue", "bitscore"]
        hit_df=pd.read_csv(hitFilepath, delimiter='\t', header=None, names=columns_lst)
        hit_df=hit_df.drop_duplicates()

        #add qfamily column
        dct_lst=[]
        for _, row in  lookup_df[["family", sstrain]].iterrows():
            if isinstance(row[sstrain],str):
                lst=row[sstrain].split(' ')
                for seqName in lst:
                    dct={}
                    dct["qfamily"]=row["family"]
                    dct["qseqid"]=seqName
                    dct_lst.append(dct)
        tmp_df=pd.DataFrame(dct_lst)
        out_df=pd.merge(hit_df, tmp_df, on="qseqid", how="left")
        assert out_df["qfamily"].isnull().sum()==0

        #add sstrain, qstrain, hit_strand columns
        dct_lst=[]
        for _,row in out_df.iterrows():
            dct={}
            dct["sstrain"]=sstrain
            dct["qstrain"]=row["qseqid"].split(':')[0]
            if row["sstart"]<row["send"]:
                dct["hit_strand"]=1
            else:
                dct["hit_strand"]=-1
            dct_lst.append(dct)
        assert out_df.shape[0]==len(dct_lst)
        tmp_df=pd.DataFrame(dct_lst, index=out_df.index)
        out_df=pd.concat([out_df, tmp_df], axis=1)
        out_df.index.name="region_id"

        outFilepath="./out/{}.csv".format(sstrain)
        out_df.to_csv(outFilepath)
        print("OUTPUT ot {}".format(outFilepath))



if __name__=="__main__":
    lookupFilepath="./out/query_lookup.csv"
    strainFilepath="/home/mitsuki/altorf/mbgd/data/ecoli/strain.lst"
    main(lookupFilepath, strainFilepath)
