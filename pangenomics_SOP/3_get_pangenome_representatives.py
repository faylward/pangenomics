import os
import sys
import subprocess
import re
import shlex
from Bio import SeqIO
from collections import defaultdict

proteinortho = open(sys.argv[1], "r")
proteins = open(sys.argv[2], "r")
ortho_out = open(sys.argv[1] +".cluster", "w")

# get protein lengths
ldict = {}
prot_dict = SeqIO.to_dict(SeqIO.parse(proteins, "fasta"))
for i in prot_dict:
	record = prot_dict[i]
	ldict[record.id] = len(record.seq)

# parse cluster IDs
tally = 0
cluster_list = []
longest_dict = defaultdict(float)
final_dict = {}
rep_list = []
for i in proteinortho.readlines():
	if i.startswith("#"):
		name = re.sub(".names.faa", "", i)
		name2 = re.sub("# ", "", name)
		ortho_out.write("Cluster\tRepresentative\tRep_Length\t"+name2)
	else:
		line = i.rstrip()
		tabs = re.split("\t|,", line)
		name = tabs[0]
		proteins = tabs[3:len(tabs)]

		max_length = defaultdict(float)
		for p in proteins:
			if p == "*":
				pass
			else:
				l = ldict[p]
				if l > max_length[p]:
					max_length[p] = l
					final_rep = p

		final_rep_length = ldict[final_rep]
		rep_list.append(final_rep)

		# get non-redundant cluster name
		tally +=1
		tally_full = str(tally).zfill(5)
		cluster = "cluster_"+ str(tally_full)
		cluster_list.append(cluster)
		ortho_out.write(cluster +"\t"+ final_rep +"\t"+ str(final_rep_length) +"\t"+ line +"\n")

record_list = []
for i in prot_dict:
	record = prot_dict[i]
	if record.id in rep_list:
		record_list.append(record)

SeqIO.write(record_list, open("representatives.faa", "w"), "fasta")
