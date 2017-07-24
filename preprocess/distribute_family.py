import pandas as pd
from Bio import SeqIO

def main(strainFilepath):
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    for strain in strain_lst:
        print("PROCESSING {}".format(strain))
        
        geneFilepath="/data/mitsuki/data/mbgd/gene/{}.gene".format(strain)
        gene_df=pd.read_csv(geneFilepath,header=None, delimiter='\t')
        
        dtype_lst=["geneseq","proteinseq"]
        for dtype in dtype_lst:
            seqFilepath="/data/mitsuki/data/mbgd/{0}/{1}.{0}".format(dtype,strain)
            name_lst=[]
            seqRec_lst=[]
            for rec in SeqIO.parse(seqFilepath, "fasta"):
                seqRec_lst.append(rec)
                name_lst.append(rec.name.split(':')[1])

            _df=pd.DataFrame(index=range(len(name_lst)))
            _df[1]=name_lst
            
            assert len(set(gene_df[1]))==gene_df.shape[0]
            df=pd.merge(_df, gene_df[[1,16]], on=1, how="left")

            for key,row in df.iterrows():
                seq=seqRec_lst[key]
                if isinstance(row[16],str):
                    family_lst=row[16].split(',')
                    for family in family_lst:
                        outFilepath="/data/mitsuki/data/mbgd/family/{0}/{1}.{0}".format(dtype,family)
                        with open(outFilepath, 'a') as f:
                            SeqIO.write(seq, f, "fasta")


if __name__=="__main__":
	strainFilepath="/home/mitsuki/altorf/mbgd/data/strain.lst"
	main(strainFilepath)
