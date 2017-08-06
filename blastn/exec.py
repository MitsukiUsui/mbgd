import pandas as pd
import numpy as np
import subprocess

def check(lookupFilepath, strainFilepath):
    lookup_df=pd.read_csv(lookupFilepath, dtype="object")
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]

    familyCount=0
    jobCount=0
    resultCount=0

    for _,row in lookup_df.iterrows():
        familyCount+=1
        msk=list(~row[strain_lst].isnull())# True if the strain does not have the family
        for i, strain in enumerate(strain_lst):
            if msk[i]:
                jobCount+=1
                filepath="/home/mitsuki/altorf/mbgd/blastn/result/{}_{}.tab".format(strain, row["family"])
                try:
                    open(filepath, 'r')
                    resultCount+=1
                except FileNotFoundError:
                    pass
    print("TOTAL")  
    print("\tfamily: {}".format(familyCount))
    print("\tjobs:   {}/{}".format(resultCount, jobCount))


def main(lookupFilepath, strainFilepath):
    lookup_df=pd.read_csv(lookupFilepath)
    strain_lst=[s.strip() for s in open(strainFilepath, 'r').readlines()]

    for _,row in lookup_df.iterrows():
        msk=list(~row[strain_lst].isnull())# True if the strain does not have the family
        for i, strain in enumerate(strain_lst):
            if msk[i]:
                cmd="/home/mitsuki/altorf/mbgd/blastn/uge_blastn_args.sh {} {} {}".format(strain,row["family"], row[strain])
                subprocess.check_call(cmd.strip().split(' '))

if __name__=="__main__":
    direc="../data/ecoli"
    strainFilepath=direc+"/strain.lst"
    lookupFilepath="query.csv"
    main(lookupFilepath, strainFilepath)
    #check(lookupFilepath, strainFilepath)
            
