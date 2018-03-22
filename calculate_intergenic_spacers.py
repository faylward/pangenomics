import os
import sys
import subprocess
import re
import shlex
from collections import defaultdict
import numpy as np

input_folder = "/home/frankaylward/Desktop/marker_gene_benchmarking/Marinimicrobia/all_genomes/gff/"
output = open("intergenic_spacers.txt", "w")

for i in os.listdir(input_folder):
	if i.endswith(".gff"):
		gff = os.path.join(input_folder, i)
		#print(i)
		infile = open(gff, "r")
		tuples = defaultdict(list)
		for j in infile.readlines():
			if j.startswith("#"):
				pass
			else:
				line = j.rstrip()
				tabs = line.split("\t")
				contig = tabs[0]
				start = int(tabs[3])
				end = int(tabs[4])

				tuples[contig].append((start, end))

		full_list = []
		for item in tuples:
			pair_list = tuples[item]
			prev_end = 0
			for n in pair_list:
				start = n[0]
				end = n[1]
				
				#if prev_end > 0:
				diff = abs(start - prev_end)
				prev_end = end

				if diff > 1000:
					pass
				else:
					full_list.append(diff)

		full_list.pop(0)
		avint = np.mean(full_list)
		output.write(i +"\t"+ str(avint) +"\n")

				


