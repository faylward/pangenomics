import os
import sys
import subprocess
import re
import shlex
import pandas
import glob
import operator
from collections import defaultdict
from Bio import SeqIO

tigr_annotations = open("/home/frankaylward/Desktop/marker_gene_benchmarking/hmm/TIGRfams/profile_annotations.tigr.txt", "r")
tigr_hmms = "/home/frankaylward/Desktop/marker_gene_benchmarking/hmm/TIGRfams/TIGRFAMs_13.0_HMM/all_tigrfam.hmm"

faa_input = sys.argv[1]
clust = sys.argv[2]
clusters = open(clust, "r")
output = open(clust+".annot", "w")

# run against all tigr
cmd = "hmmsearch --cpu 16 --cut_nc --tblout representative.hmmout "+ tigr_hmms +" "+ faa_input
cmd2 = shlex.split(cmd)
#subprocess.call(cmd2, stdout=open("hmm.out", 'a'), stderr=open("error_file.txt", 'a'))

# loop through the hmm output and get the best hits
f = open("representative.hmmout", 'r')
hit_dict = defaultdict(lambda:"no_hit")
bit_dict = defaultdict(lambda:float(0))
hit_type = {}
marker_dict = {}

for line in f.readlines():
	if line.startswith("#"):
		pass
	else:
		newline = re.sub( '\s+', '\t', line)
		list1 = newline.split('\t')
		ids = list1[0]
		hit = list1[2]
		acc2 = list1[3]
		bit_score = list1[5]
		score = float(bit_score)
				
		if score > bit_dict[ids]:
			hit_dict[ids] = hit
			bit_dict[ids] = score

# get an annotation dictionary
annot_dict = defaultdict(lambda:"no_hit")
for i in tigr_annotations:
	line = i.rstrip()
	tabs = line.split("\t")
	tigr = tabs[0]
	annot = tabs[3]
	annot_dict[tigr] = annot

# compile the final table
for i in clusters.readlines():
	if i.startswith("Cluster\t"):
		line = i.rstrip()
		tabs = line.split("\t")
		print tabs[0] +"\trep\tannot\tscore\thit_desc\t"+ "\t".join(tabs[2:len(tabs)]) +"\n"
		output.write(tabs[0] +"\trep\tannot\tscore\t"+ "\t".join(tabs[2:len(tabs)]) +"\n")
	
	else:
		line = i.rstrip()
		tabs = line.split("\t")
		rep = tabs[1]
		annot = hit_dict[rep]
		score = str(bit_dict[rep])	
		output.write(tabs[0] +"\t"+ rep +"\t"+ annot +"\t"+ score +"\t"+ annot_dict[annot] +"\t"+  "\t".join(tabs[2:len(tabs)]) +"\n")







