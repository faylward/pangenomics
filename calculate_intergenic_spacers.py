import os
import sys
import subprocess
import re
import shlex
from collections import defaultdict
import numpy as np

input_folder = sys.argv[1]
output = open(sys.argv[2], "w")
output.write("genome\tnum_spacers\tmean_len\tmedian_len\tminimum_len\tmaximum_len\n")

for i in os.listdir(input_folder):
	if i.endswith(".gff"):
		gff = os.path.join(input_folder, i)
		print(i)
		infile = open(gff, "r")
		tuples = defaultdict(list)
		for j in infile.readlines():
			if j.startswith("##FASTA"):
				break
			elif j.startswith("#"):
				pass
			else:
				line = j.rstrip()
				tabs = line.split("\t")
				#print(tabs)
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
				
				diff = start - prev_end
				
				if diff < 0:
					full_list.append(0)
				else:
					full_list.append(diff)
					
				prev_end = end
				#full_list.append(diff)

		full_list.pop(0)
		avint = np.mean(full_list)
		median = np.median(full_list)
		mini = min(full_list)
		maxi = max(full_list)
		#print(full_list)
		output.write(i +"\t"+ str(len(full_list)) +"\t"+ str(round(avint, 1)) +"\t"+ str(median) +"\t"+ str(mini) +"\t"+ str(maxi) +"\n")

