import pandas as pd

def main(strainFilepath):
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    for strain in strain_lst:
        filepath="/data/mitsuki/data/mbgd/gene/{}.gene".format(strain)
        gene_df=pd.read_csv(filepath, delimiter='\t', header=None)
        
        out_df=pd.DataFrame(index=gene_df.index, columns=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
        out_df["seqname"]=gene_df[0]+':'+gene_df[15]
        out_df["source"]="MBGD"
        out_df["feature"]=gene_df[7]
        out_df["start"]=gene_df[3]
        out_df["end"]=gene_df[4]
        out_df["score"]="."
        out_df.loc[gene_df[5]==1, "strand"]="+"
        out_df.loc[gene_df[5]==-1, "strand"]="-"
        out_df["frame"]="."
        out_df["attribute"]="name "+gene_df[16]+";gene_name "+gene_df[1].fillna('')+';description '+gene_df[13].fillna('')

        outFilepath="/data/mitsuki/data/mbgd/gff/{}.gff".format(strain)
        out_df.to_csv(outFilepath, index=False, header=None, sep='\t')
    
if __name__=="__main__":
    strainFilepath="/home/mitsuki/altorf/mbgd/data/strain.lst"
    main(strainFilepath)
