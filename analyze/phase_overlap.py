#!/home/mitsuki/.pyenv/versions/anaconda3-4.3.1/bin/python

import sys
import pandas as pd
import math
import numpy as np
from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo as matlist

from phase import Phase

matrix = matlist.blosum62


def main(seqFilepath, overlapFilepath, outFilepath):
    seqRec_dct={}
    for seqRec in SeqIO.parse(seqFilepath, "fasta"):
        seqRec_dct[seqRec.id]=seqRec
    overlap_df=pd.read_csv(overlapFilepath)

    dct_lst=[]
    for overlapId in overlap_df["overlap_id"]:
        print(overlapId)
        dct={}
        dct["overlap_id"]=overlapId
        
        qseq_dna=seqRec_dct["{}:qseq_dna".format(overlapId)].seq
        sseq_dna=seqRec_dct["{}:sseq_dna".format(overlapId)].seq
        qseq_pro=seqRec_dct["{}:qseq_pro".format(overlapId)].seq
        sseq_pro=seqRec_dct["{}:sseq_pro".format(overlapId)].seq
       
        # get 6 frame translation of qseq_dna
        qtrans_lst=[]
        for strand in [1,-1]:
            for gap in range(3):
                start=gap
                end=gap+int(math.floor(len(qseq_dna)-gap)/3)*3

                seq=qseq_dna[gap:end]
                if strand==-1:
                    seq=seq.reverse_complement()
                qtrans_lst.append(str(seq.translate(table=11)).replace('*','X'))
        assert len(qtrans_lst)==6
            
        qscore_lst=[]
        sscore_lst=[]
        for qtrans in qtrans_lst:
            alns_pro=pairwise2.align.globalds(qtrans, qseq_pro, matrix, -10, -0.5)
            qscore_lst.append(alns_pro[0][2])

            alns_pro=pairwise2.align.globalds(qtrans, sseq_pro, matrix, -10, -0.5)
            sscore_lst.append(alns_pro[0][2])
            
        for i in range(6):
            dct["qscore{}".format(i)]=qscore_lst[i]
            dct["sscore{}".format(i)]=sscore_lst[i]
        phase=Phase()
        dct["relative"]=phase.phase_lst[phase.relative_int(np.argmax(qscore_lst), np.argmax(sscore_lst))]
        dct_lst.append(dct)
    
    out_df=pd.DataFrame(dct_lst)
    columns_lst=["overlap_id"]
    for i in range(6):
        columns_lst.append("qscore{}".format(i))
    for i in range(6):
        columns_lst.append("sscore{}".format(i))
    columns_lst.append("relative")
    out_df=out_df[columns_lst]
    out_df.to_csv(outFilepath, index=False)
    print("OUTPUT to {}".format(outFilepath))


if __name__=="__main__":
    #strain="eab"
    strain=sys.argv[1]
    direc="/home/mitsuki/altorf/mbgd/analyze/out"
    seqFilepath=direc+"/{}_ovr.fasta".format(strain)
    overlapFilepath=direc+"/{}_filtered.csv".format(strain)
    outFilepath=direc+"/{}_relative.csv".format(strain)
    main(seqFilepath, overlapFilepath, outFilepath)
