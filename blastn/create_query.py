#!/usr/bin/env python3

import subprocess
import os
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extract_sequence(orfId, dtype="geneseq", family_tpl=None, start=None, end=None):
    """
    argument:
        orfId...ecg:E2348C_1809
        dtype...geneseq or proteinseq
        family_tpl...(sulfolobus_islandicus ,family1) as e.g.
    """
    
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

def main(target, subStrain, lookupFilepath):
    lookup_df=pd.read_csv(lookupFilepath, index_col=0)
    
    print("START: extract query sequences for {}".format(subStrain))
    seqRec_lst=[]
    for family, orfId_lst in lookup_df[subStrain].dropna().iteritems():
        for orfId in orfId_lst.split(" "):
            seq=extract_sequence(orfId, family_tpl=(target, family))
            seqRec=SeqRecord(seq, id=orfId, description="")
            seqRec_lst.append(seqRec)
        
    queryDirec="./query/{}".format(target)
    os.makedirs(queryDirec, exist_ok=True)
    queryFilepath="{}/{}.query".format(queryDirec, subStrain)
    with open(queryFilepath, 'w') as f:
        for seqRec in seqRec_lst:
            SeqIO.write(seqRec, f, "fasta")
    print("DONE: output {} querys to {}".format(len(seqRec_lst), queryFilepath))


if __name__=="__main__":
    target=sys.argv[1]
    subStrain=sys.argv[2]
    lookupFilepath="../data/{}/query_lookup.csv".format(target)
    main(target, subStrain, lookupFilepath)
