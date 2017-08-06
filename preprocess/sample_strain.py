import random
random.seed(42)

def main(clusterFilepath, outFilepath, n=None):
    with open(clusterFilepath, 'r') as f:
        column_lst=f.readline().strip().split('\t')
    all_lst=column_lst[3:-8]
    if n is None:
        strain_lst=all_lst
    else:
        strain_lst=random.sample(all_lst, n)
    
    with open(outFilepath, 'w') as f:
        for strain in strain_lst:
            f.write(strain+'\n')
    
    print("SAMPLED {}/{} to {}".format(len(strain_lst), len(all_lst), outFilepath))
    
if __name__=="__main__":
    direc="../data/ecoli"
    clusterFilepath=direc+"/ecoli_cluster.tab"
    outFilepath=direc+"/strain.lst"
    main(clusterFilepath, outFilepath)
    
