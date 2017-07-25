#!/home/mitsuki/.pyenv/versions/anaconda3-4.3.1/bin/python

import pandas as pd
import numpy as np
import sys
import subprocess
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SubsMat import MatrixInfo as matlist

def extract_sequence(orfId, start=None, end=None, dtype="geneseq"):
	"""
	argument:
		orfId...ecg:E2348C_1809
		dtype...geneseq or proteinseq
	"""
	
	dataDir="/data/mitsuki/data/mbgd/"+dtype
	cmd="/home/mitsuki/usr/bin/fatt extract --seq {0} {1}/{2}.{3}".format(orfId, dataDir, orfId.split(':')[0], dtype)
	result=subprocess.check_output(cmd.strip().split(' ')).decode('utf-8')
	result=''.join(result.split('\n')[1:])
	
	if result=='':
		print("NOT FOUND: {}".format(orfId))
		return None
	else:
		if start is None:
			start=0
		else:
			start=max(0,start)
		
		if end is None :
			if dtype=="proteinseq":
				end=len(result)-1
			else:
				end=len(result)
		else:
			if dtype=="proteinseq":
				end=min(len(result)-1,end)
			else:
				end=min(len(result),end)
		return Seq(result[start:end])

def main(strain):
	strainFilepath =  "/home/mitsuki/altorf/mbgd/blastn/out/{}.csv".format(strain)
	overlapFilepath = "/home/mitsuki/altorf/mbgd/analyze/out/{}_ovr.csv".format(strain)
	logFilepath =     "/home/mitsuki/altorf/mbgd/analyze/out/{}_log.csv".format(strain)
	outFilepath =     "/home/mitsuki/altorf/mbgd/analyze/out/{}_out.csv".format(strain)

	strain_df=pd.read_csv(strainFilepath)
	overlap_df=pd.read_csv(overlapFilepath)
	merged_df=pd.merge(overlap_df, strain_df, how="left", on="region_id")
	
	merged_df["hit_strand"]=1
	merged_df.loc[(merged_df["sstart"]>merged_df["send"]), "hit_strand"]=-1

	###
	# 1. calcurate the relative position of overlapping region to blastn hit, denoted by start_per, end_per each
	###
	tpl_lst=[]
	for _,row in merged_df.iterrows():
		length=abs(row["send"]-row["sstart"])
		if row["hit_strand"]==1:
			start_per=(row["ofirst"]-row["sstart"])/length
			end_per=(row["olast"]-row["sstart"])/length
		else:
			start_per=(row["sstart"]-row["olast"])/length
			end_per=(row["sstart"]-row["ofirst"])/length
		tpl_lst.append((start_per,end_per))
	merged_df["start_per"]=0
	merged_df["end_per"]=0
	merged_df[["start_per","end_per"]]=tpl_lst
	
	###
	# 2. calculate qfirst_pro, qlast_pro
	###
	merged_df["qstart_gen"]=np.ceil(merged_df["qstart"]+(merged_df["qend"]-merged_df["qstart"])*merged_df["start_per"]).astype(int)
	merged_df["qend_gen"] =np.floor(merged_df["qstart"]+(merged_df["qend"]-merged_df["qstart"])*merged_df["end_per"]).astype(int)
	merged_df["qstart_pro"]=np.ceil(merged_df["qstart_gen"]/3).astype(int)
	merged_df["qend_pro"]=np.floor(merged_df["qend_gen"]/3).astype(int)
	

	###
	# 3. calc alignment score for geneseq & proteinseq
	###
	matrix = matlist.blosum62
	
	scoreGen_lst=[]
	scorePro_lst=[]

	print("START {} alignments".format(merged_df.shape[0]))
	with open(logFilepath, 'w') as f:
		for key,row in merged_df.iterrows():
			if key%100==0:
				print("\t{}".format(key))
			
			if row["olength"]<10: #to short for alignment
				scoreGen_lst.append(0)
				scorePro_lst.append(0)
			elif row["olength"]>=10:
				f.write("***{}***\n".format(key))
				
				#geneseq alignment
				qseq_gen=extract_sequence(row["qseqid"], 
										  row["qstart_gen"], row["qend_gen"], dtype="geneseq")
				sseq_gen=extract_sequence(row["chr_name"].split(':')[0]+':'+row["cds_name"],
										  row["sstart_gen"], row["send_gen"], dtype="geneseq")
				if row["hit_strand"]*row["cds_strand"]==1:
					alns_gen=pairwise2.align.globalms(qseq_gen, sseq_gen, 2, -1, -.5, -.1)
				elif row["hit_strand"]*row["cds_strand"]==-1:
					alns_gen=pairwise2.align.globalms(qseq_gen, sseq_gen.reverse_complement(), 2, -1, -.5, -.1)
				scoreGen_lst.append(alns_gen[0][2])

				f.write(str(qseq_gen)+'\n')
				f.write(str(sseq_gen)+'\n')
				f.write(str(scoreGen_lst[-1])+'\n')
			

				#proteinseq alignment
				qseq_pro=extract_sequence(row["qseqid"], 
										  row["qstart_pro"], row["qend_pro"], dtype="proteinseq")
				sseq_pro=extract_sequence(row["chr_name"].split(':')[0]+':'+row["cds_name"],
										  row["sstart_pro"], row["send_pro"], dtype="proteinseq")
				alns_pro=pairwise2.align.globalds(qseq_pro, sseq_pro, matrix, -10, -0.5)
				scorePro_lst.append(alns_pro[0][2])
				
				f.write(str(qseq_pro)+'\n')
				f.write(str(sseq_pro)+'\n')
				f.write(str(scorePro_lst[-1])+'\n\n')
	merged_df["score_gen"]=scoreGen_lst
	merged_df["score_pro"]=scorePro_lst
	merged_df.to_csv(outFilepath,index=False)

if __name__=="__main__":
	strain=sys.argv[1]
	main(strain)
	
