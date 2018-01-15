import pandas as pd

def main(strainFilepath, outFilepath):
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    strain_lst.remove("gm04124")
    strain_lst.remove("gm04125")
    all_df=pd.DataFrame(columns=["overlap_id", "qfamily", "sfamily", "relative"])

    for strain in strain_lst:
        print(strain)
        overlapFilepath="./out/{}_filtered.csv".format(strain)
        relativeFilepath="./out/{}_relative.csv".format(strain)

        overlap_df=pd.read_csv(overlapFilepath)
        relative_df=pd.read_csv(relativeFilepath, dtype={"relative":str})
    
        #break sfamily_lst to each row
        dct_lst=[]
        for _, row in overlap_df.iterrows():
            if isinstance(row["sfamily"],str):
                for sfamily in row["sfamily"].split(','):
                    dct={}
                    dct["overlap_id"]=row["overlap_id"]
                    dct["qfamily"]=row["qfamily"]
                    dct["sfamily"]=sfamily
                    dct_lst.append(dct)
        tmp_df=pd.DataFrame(dct_lst)
        print(tmp_df.shape)
        tmp_df.head()
        
        dot_df=pd.merge(tmp_df, relative_df[["overlap_id", "relative"]])
        assert dot_df.shape[0]==tmp_df.shape[0]
        all_df=pd.concat([all_df, dot_df])
    
    all_df=all_df[["qfamily", "sfamily", "relative"]] 
    out_df = pd.DataFrame(all_df.groupby(["qfamily", "sfamily", "relative"]).size().rename('weight'))
    out_df.to_csv(outFilepath)

if __name__=="__main__":
    strainFilepath="/home/mitsuki/altorf/mbgd/data/ecoli/strain.lst"
    outFilepath="./out/phase_dot.csv"
    main(strainFilepath, outFilepath)
