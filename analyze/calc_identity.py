#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import subprocess
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SubsMat import MatrixInfo as matlist

def extract_sequence(orfId, dtype="geneseq", family_tpl=None, start=None, end=None, ):
    """
    argument:
        orfId...ecg:E2348C_1809
        dtype...geneseq or proteinseq
        family_tpl...(target, familyName). ("sulfolobus_islandicus" ,"family1") as e.g.
    """
    assert dtype in ("geneseq", "proteinseq")
    
    if family_tpl is None:
        seqFilepath="/data/mitsuki/data/mbgd/{0}/{1}.{0}".format(dtype, orfId.split(':')[0])
    else:
        strain=orfId.split(':')[0]
        seqFilepath="/data/mitsuki/data/mbgd/{0}/family/{1}/{2}.{1}".format(family_tpl[0], dtype, family_tpl[1])
        
    cmd="fatt extract --seq {0} {1}".format(orfId, seqFilepath)
    result=subprocess.check_output(cmd.strip().split(' ')).decode('utf-8')
    result=''.join(result.split('\n')[1:])
    
    if result=='':
        print("NOT FOUND: {}".format(orfId))
        print("\t{}".format(seqFilepath))
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

def main(target, overlapFilepath, logFilepath):
    overlap_df=pd.read_csv(overlapFilepath, index_col=0)
    
    # calc alignment score for geneseq & proteinseq
    matrix = matlist.blosum62
    scoreDna_lst=[]
    scorePro_lst=[]
    print("TOTAL {} alignments".format(overlap_df.shape[0]), flush = True)
    with open(logFilepath, 'w') as f:
        for key,row in overlap_df.iterrows():
            if key%100==0:
                print("\t{}".format(key), flush = True)

            if row["olength"]<10: #to short for alignment
                scoreDna_lst.append(0)
                scorePro_lst.append(0)
            elif row["olength"]>=10:
                #geneseq alignment
                qseq_dna=extract_sequence(row["qorf_id"],  dtype="geneseq", family_tpl=(target, row["qfamily"]),
                                                    start=row["qstart_dna"], end = row["qend_dna"])
                sseq_dna=extract_sequence(row["sorf_id"],  dtype="geneseq",
                                                    start=row["sstart_dna"], end = row["send_dna"])
                
                try:
                    if row["hit_strand"]*row["cds_strand"]==1:
                        alns_dna=pairwise2.align.globalms(qseq_dna, sseq_dna, 2, -1, -.5, -.1)
                    elif row["hit_strand"]*row["cds_strand"]==-1:
                        alns_dna=pairwise2.align.globalms(qseq_dna, sseq_dna.reverse_complement(), 2, -1, -.5, -.1)
                    scoreDna_lst.append(alns_dna[0][2])
                except (IndexError, SystemError):
                    scoreDna_lst.append(0)

                f.write(">{}:{}\n".format(key, "qseq_dna"))
                f.write(str(qseq_dna)+'\n')
                f.write(">{}:{}\n".format(key, "sseq_dna"))
                f.write(str(sseq_dna)+'\n')
            
                #proteinseq alignment
                qseq_pro=extract_sequence(row["qorf_id"],  dtype="proteinseq", family_tpl=(target, row["qfamily"]),
                                                    start=row["qstart_pro"], end = row["qend_pro"])
                sseq_pro=extract_sequence(row["sorf_id"],  dtype="proteinseq",
                                                    start=row["sstart_pro"], end = row["send_pro"])
                qseq_pro=str(qseq_pro).replace('*','X')
                sseq_pro=str(sseq_pro).replace('*','X')
                
                try: 
                    alns_pro=pairwise2.align.globalds(qseq_pro, sseq_pro, matrix, -10, -0.5)
                    scorePro_lst.append(alns_pro[0][2])
                except (IndexError, SystemError):
                    scorePro_lst.append(0)

                f.write(">{}:{}\n".format(key, "qseq_pro"))
                f.write(str(qseq_pro)+'\n')
                f.write(">{}:{}\n".format(key, "sseq_pro"))
                f.write(str(sseq_pro)+'\n')
                    
    overlap_df["score_dna"]=scoreDna_lst
    overlap_df["score_pro"]=scorePro_lst
    overlap_df.to_csv(overlapFilepath)

if __name__=="__main__":
    target=sys.argv[1]
    strain=sys.argv[2]
    overlapFilepath="./out/{}/{}_ovr.csv".format(target, strain)
    logFilepath="./out/{}/{}_ovr.fasta".format(target, strain)
    main(target, overlapFilepath, logFilepath)
    
