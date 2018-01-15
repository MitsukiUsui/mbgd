#!/usr/bin/env python3

import sys

target=sys.argv[1]
direc="/home/mitsuki/altorf/mbgd/data/{}".format(target)
clusterFilepath="{}/cluster.tab".format(direc)
strainFilepath="{}/strain.lst".format(direc)
with open(clusterFilepath, 'r') as f:
    column_lst=f.readline().strip().split('\t')
strain_lst=column_lst[3:-6]

with open(strainFilepath, "w") as f:
    for strain in strain_lst:
        f.write("{}\n".format(strain))
print("DONE: output {} strains to {}".format(len(strain_lst), strainFilepath))