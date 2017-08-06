#!/home/mitsuki/.pyenv/versions/anaconda3-4.3.1/bin/python

import pandas as pd
import numpy as np
import sys
import subprocess
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SubsMat import MatrixInfo as matlist


def extract_sequence(orfId, start=None, end=None, dtype="geneseq", family=None):
    """
    argument:
        orfId...ecg:E2348C_1809
        dtype...geneseq or proteinseq
    """
    
    if family is None:
        seqFilepath="/data/mitsuki/data/mbgd/{0}/{1}.{0}".format(dtype, orfId.split(':')[0])
    else:
        strain=orfId.split(':')[0]
        seqFilepath="/data/mitsuki/data/mbgd/family/{0}/{1}/{1}_{2}.{0}".format(dtype, family, strain)

    cmd="/home/mitsuki/usr/bin/fatt extract --seq {0} {1}".format(orfId, seqFilepath)
    result=subprocess.check_output(cmd.strip().split(' ')).decode('utf-8')
    result=''.join(result.split('\n')[1:])
    
    if result=='':
        print("NOT FOUND: {}".format(orfId))
        print(seqFilepath)
        return None
    else:
        if start is None:
            start=0
        else:
            start=max(0,start)
        
        if end is None:
            end=len(result)
        else:
            end=min(len(result),end)
        return Seq(result[start:end])

def main(overlapFilepath, logFilepath):
    overlap_df=pd.read_csv(overlapFilepath, index_col=0)
    
    # calc alignment score for geneseq & proteinseq
    matrix = matlist.blosum62
    scoreDna_lst=[]
    scorePro_lst=[]
    print("TOTAL {} alignments".format(overlap_df.shape[0]))
    with open(logFilepath, 'w') as f:
        for key,row in overlap_df.iterrows():
            if key%100==0:
                print("\t{}".format(key))

            if row["olength"]<10: #to short for alignment
                scoreDna_lst.append(0)
                scorePro_lst.append(0)
            elif row["olength"]>=10:
                
                #geneseq alignment
                qseq_dna=extract_sequence(row["qorf_id"], 
                                          row["qstart_dna"], row["qend_dna"], dtype="geneseq", family=row["qfamily"])
                sseq_dna=extract_sequence(row["sorf_id"],
                                          row["sstart_dna"], row["send_dna"], dtype="geneseq")
                if row["hit_strand"]*row["cds_strand"]==1:
                    alns_dna=pairwise2.align.globalms(qseq_dna, sseq_dna, 2, -1, -.5, -.1)
                elif row["hit_strand"]*row["cds_strand"]==-1:
                    alns_dna=pairwise2.align.globalms(qseq_dna, sseq_dna.reverse_complement(), 2, -1, -.5, -.1)
                scoreDna_lst.append(alns_dna[0][2])

                f.write(">{}:{}\n".format(key, "qseq_dna"))
                f.write(str(qseq_dna)+'\n')
                f.write(">{}:{}\n".format(key, "sseq_dna"))
                f.write(str(sseq_dna)+'\n')
            
                #proteinseq alignment
                qseq_pro=extract_sequence(row["qorf_id"], 
                                          row["qstart_pro"], row["qend_pro"], dtype="proteinseq", family=row["qfamily"])
                sseq_pro=extract_sequence(row["sorf_id"],
                                          row["sstart_pro"], row["send_pro"], dtype="proteinseq")
                qseq_pro=str(qseq_pro).replace('*','X')
                sseq_pro=str(sseq_pro).replace('*','X')
                alns_pro=pairwise2.align.globalds(qseq_pro, sseq_pro, matrix, -10, -0.5)
                scorePro_lst.append(alns_pro[0][2])

                f.write(">{}:{}\n".format(key, "qseq_pro"))
                f.write(str(qseq_pro)+'\n')
                f.write(">{}:{}\n".format(key, "sseq_pro"))
                f.write(str(sseq_pro)+'\n')
    overlap_df["score_dna"]=scoreDna_lst
    overlap_df["score_pro"]=scorePro_lst
    overlap_df.to_csv(overlapFilepath)

if __name__=="__main__":
    strain=sys.argv[1]
    overlapFilepath = "/home/mitsuki/altorf/mbgd/analyze/out/{}_ovr.csv".format(strain)
    logFilepath =     "/home/mitsuki/altorf/mbgd/analyze/out/{}_ovr.fasta".format(strain)
    main(overlapFilepath, logFilepath)
    
