import pandas as pd
import subprocess
import re

def extract_sequence(orfId, dtype):
    """
    argument:
        orfId...ecg:E2348C_1809
        dtyp...geneseq or proteinseq
    """
    
    dataDir="/data/mitsuki/data/mbgd/"+dtype
    cmd="fatt extract --seq {0} {1}/{2}.{3}".format(orfId, dataDir, orfId.split(':')[0], dtype)
    result=subprocess.check_output(cmd.strip().split(' ')).decode('utf-8')
    result=''.join(result.split('\n')[1:])
    
    if result=='':
        print("NOT FOUND: {}".format(orfId))
    
    return result


def main(strainFilepath):
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]
    df=pd.read_csv("../preprocess/sampled_cluster.csv")
    pattern = r"([^()]+)(\([0-9]+\))?"
    r=re.compile(pattern)

    for key, row in df.iterrows():
        if key%10==0:
            print("PROCESSING {}".format(key))

        outFilepath="./query/{}.fna".format(row["Family"])
        with open(outFilepath, 'w') as f:
            for k,v in row[strain_lst].dropna().iteritems():
                for orfId in v.split():
                    orfId=r.findall(orfId)[0][0]
                    header='>'+orfId
                    f.write(header+'\n')
                    f.write(extract_sequence(orfId, dtype="geneseq")+'\n')


if __name__=="__main__":
    main("../preprocess/strain.lst")