import pandas as pd

def main(strainFilepath):
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    for strain in strain_lst:
        overlapFilepath="/home/mitsuki/altorf/mbgd/analyze/out/{}_ovr.csv".format(strain)
        outFilepath="/home/mitsuki/altorf/mbgd/analyze/out/{}_filtered.csv".format(strain)
        overlap_df=pd.read_csv(overlapFilepath)
        out_df=overlap_df[(overlap_df["olength"]>200) & (overlap_df["score_pro"]<0)]
        out_df.to_csv(outFilepath, index=False)
        print("OUTPUT to {}".format(outFilepath))

if __name__=="__main__":
    strainFilepath="/home/mitsuki/altorf/mbgd/data/ecoli/strain.lst"
    #strainFilepath="strain.lst"
    main(strainFilepath)
