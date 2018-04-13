import os
import sys
import subprocess
import re
import shlex
from Bio import SeqIO
from collections import defaultdict

finput = open("/home/frankaylward/Desktop/marker_gene_benchmarking/Marinimicrobia/functional_annotations/final_set_dir.proteinortho", "r")
lengths = open("/home/frankaylward/Desktop/marker_gene_benchmarking/Marinimicrobia/functional_annotations/protein_lengths.txt", "r")
longest_rep = open("/home/frankaylward/Desktop/marker_gene_benchmarking/Marinimicrobia/functional_annotations/longest_rep.txt", "w")

pangenome = open("/home/frankaylward/Desktop/marker_gene_benchmarking/Marinimicrobia/functional_annotations/pangenome.faa", "r")
ortho_out = open("/home/frankaylward/Desktop/marker_gene_benchmarking/Marinimicrobia/functional_annotations/final_set_dir.proteinortho.cluster", "w")

ldict = {}
for i in lengths.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	name = tabs[0]
	length = tabs[1]
	ldict[name] = float(length)

tally = 0
cluster_list = []
longest_dict = defaultdict(float)
final_dict = {}
for i in finput.readlines():
	if i.startswith("#"):
		name = re.sub("_genomic.prodigal.faa", "", i)
		name1 = re.sub("_protein.faa", "", name)
		name2 = re.sub(".faa", "", name1)
		name3 = re.sub("\# ", "", name2)
		ortho_out.write("Cluster\t"+name3)
	else:
		line = i.rstrip()
		#tabs = line.split("\t")
		tabs = re.split("\t|,", line)
		name = tabs[0]
		proteins = tabs[3:len(tabs)]

		tally +=1
		tally_full = str(tally).zfill(5)
		cluster = "cluster_"+ str(tally_full)
		cluster_list.append(cluster)

		ortho_out.write(cluster +"\t"+ line +"\n")

		#print line
		#print proteins

		for j in proteins:
			if j == "*":
				pass
			else:
				#print cluster, j, ldict[j]
				if ldict[j] > longest_dict[cluster]:
					final_dict[cluster] = j
					longest_dict[cluster] = ldict[j]

final_protein_set = []
for i in cluster_list:
	final_protein_set.append(final_dict[i])
	longest_rep.write(i +"\t"+ final_dict[i] +"\t"+ str(ldict[final_dict[i]]) +"\n")


for record in SeqIO.parse(pangenome, "fasta"):
	if record.id in final_protein_set:
		#print record.id, "all good"
		pass
	else:
		print record.id

