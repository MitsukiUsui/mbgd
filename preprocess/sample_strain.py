import random
random.seed(42)

def main(clusterFilepath, outFilepath, n=10):
    with open(clusterFilepath, 'r') as f:
        column_lst=f.readline().strip().split('\t')
    all_lst=column_lst[3:-8]
    strain_lst=random.sample(all_lst, n)
    
    with open(outFilepath, 'w') as f:
        for strain in strain_lst:
            f.write(strain+'\n')
    
    print("SAMPLED {}/{} to {}".format(n, len(all_lst), outFilepath))
    
if __name__=="__main__":
    direc="../data/streptomyces"
    clusterFilepath=direc+"/streptomyces_cluster.tab"
    outFilepath=direc+"/strain.lst"
    n=10
    main(clusterFilepath, outFilepath, n)
    
