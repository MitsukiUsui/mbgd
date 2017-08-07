import pandas as pd
from collections import defaultdict
from graphviz import Digraph
from phase import Phase

def get_comp_ddct(dot_df):
    """
    decompose to connected component using union-find
    """

    def get_root(root_lst, i):
    if root_lst[i]==i:
        return root_lst, i
    else:
        root_lst, root=get_root(root_lst, root_lst[i])
        root_lst[i]=root
        return root_lst, root
    
    #get node_lst
    node_set=set()
    node_set.update(dot_df["qfamily"])
    node_set.update(dot_df["sfamily"])
    node_lst=list(node_set)
    
    root_lst=list(range(len(node_lst)))#initialize
    for _,row in dot_df.iterrows(): #process each edge
        qidx=node_lst.index(row["qfamily"])
        sidx=node_lst.index(row["sfamily"])
        root_lst, qroot=get_root(root_lst, qidx)
        root_lst, proot=get_root(root_lst, sidx)
        root_lst[qroot]=min(qroot, proot)
        root_lst[proot]=min(qroot, proot)
    for i in range(len(node_lst)):#finalization process to point each root
        root_lst, _=get_root(root_lst, i)
    
    #distibute node based on root
    comp_ddct=defaultdict(list)
    for i in range(len(node_lst)):
        comp_ddct[root_lst[i]].append(node_lst[i])
        
    return comp_ddct


def get_event_df(dot_df, comp_ddct):
    phase=Phase()
    ddct_lst=[]
    for _, family_lst in comp_ddct.items():
        assert len(family_lst)>=2
        tmp_df=dot_df[dot_df["qfamily"].isin(family_lst)]
        f2p={} #family to phase

        for _, row in tmp_df.iterrows():
            if (row["qfamily"] in f2p.keys()) and (row["sfamily"] in f2p.keys()):
                qp=f2p[row["qfamily"]]
                sp=f2p[row["sfamily"]]
                assert phase.relative(qp, sp)==row["relative"]
            elif (row["qfamily"] in f2p.keys()) and not(row["sfamily"] in f2p.keys()):
                qp=f2p[row["qfamily"]]
                sp=phase.operate(qp, row["relative"])
                f2p[row["sfamily"]]=sp
            elif not(row["qfamily"] in f2p.keys()) and (row["sfamily"] in f2p.keys()):
                sp=f2p[row["sfamily"]]
                revop=phase.relative(row["relative"], "+0")# revops の関数をphaseに作るべき
                qp=phase.operate(sp, revop)
                f2p[row["qfamily"]]=qp
            else:
                f2p[row["qfamily"]]="+0"
                f2p[row["sfamily"]]=row["relative"]

        ddct=defaultdict(list)
        for f, p in f2p.items():
            ddct[p].append(f)
        ddct_lst.append(ddct)

    dct_lst=[]
    for ddct in ddct_lst:
        dct={}
        for k,v in ddct.items():
            dct[k]=",".join(v)
        dct_lst.append(dct)
    out_df=pd.DataFrame(dct_lst)
    return out_df

def main(overlapFilepath, relativeFilepath):
    overlap_df=pd.read_csv(overlapFilepath)
    relative_df=pd.read_csv(relativeFilepath, dtype={"relative":str})
    
    #break sfamily_lst to each row
    dct_lst=[]
    for _, row in filtered_df.iterrows():
        for sfamily in row["sfamily"].split(','):
            dct={}
            dct["overlap_id"]=row["overlap_id"]
            dct["qfamily"]=row["qfamily"]
            dct["sfamily"]=sfamily
            dct_lst.append(dct)
    tmp_df=pd.DataFrame(dct_lst)
    print(tmp_df.shape)
    tmp_df.head()
    
    dot_df=pd.merge(tmp_df, relative_df[["overlap_id", "relative"]])
    assert dot_df.shape[0]==tmp_df.shape[0]
    comp_ddct=get_comp_ddct(dot_df)
    out_df=get_event_df(dot_df, comp_ddct)
    

if __name__=="__main__":
    overlapFilepath="./test/eab_filtered.csv"
    relativeFilepath="./test/eab_relative.csv"
    main(overlapFilepath, relativeFilepath)