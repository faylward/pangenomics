import os
import sys
import re
from collections import defaultdict
import operator
import subprocess
import shlex
import numpy as np

proteinortho = open("test_genome.proteinortho", "r")
output = open("genome_distances.txt", "w")
output.write("genome_pair\tnum_prot1\tnum_prot2\tnum_shared\tprop1\tprop2\tav_prop\n")

cluster_dict = defaultdict(list)
tally = int(0)
for i in proteinortho.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	entries = tabs[3:len(tabs)]

	if line.startswith("#"):
		genomes = [re.sub("_protein.faa", "", n) for n in entries]

	else:
		tally +=1
		cluster = "cluster_"+str(tally)

		for index, protein in enumerate(entries):
			#print index, protein
			if protein == "*":
				pass
			else:
				genome = genomes[index]
				cluster_dict[genome].append(cluster)

already_done = {}
for i in cluster_dict:
	for j in cluster_dict:
		comparison = i +"||"+ j
		inverse = j +"||"+ i

		if inverse in already_done:
			pass
		else:
			already_done[comparison] = comparison

			genome1 = cluster_dict[i]
			genome2 = cluster_dict[j]
			g1length = len(genome1)
			g2length = len(genome2)
		
			shared_proteins = list(set(genome1).intersection(genome2))
			num_shared = len(shared_proteins)

			prop1 = float(num_shared) / float(g1length)
			prop2 = float(num_shared) / float(g2length)
			mean = np.mean([prop1, prop2])

			output.write(comparison +"\t"+ str(g1length) +"\t"+ str(g2length) +"\t"+ str(num_shared) +"\t"+ str(prop1) +"\t"+ str(prop2) +"\t"+ str(mean) +"\n")
