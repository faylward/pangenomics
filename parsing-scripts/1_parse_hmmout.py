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

finput = "/home/frankaylward/Desktop/marker_gene_benchmarking/Marinimicrobia/functional_annotations/hmmer_annotations"

################################################################
######### Loop through and parse the HMM outputs ###############
################################################################

for i in os.listdir(finput):

	if i.endswith(".tblout"):

		print i
		hits = []
		bit_dict = {}
		len_dict = {}
		hit_dict = {}
		bit_dict = defaultdict(int)
		hit_type = {}

		hmm_in = os.path.join(finput, i)
		hmmer = open(hmm_in, "r")
		out = re.sub(".tblout", ".parsed", hmm_in)
		output = open(out, "w")

		for j in hmmer.readlines():
			if j.startswith("#"):
				pass
			else:
				line = j.rstrip()
				newline = re.sub( '\s+', '\t', line)
				list1 = newline.split('\t')
				query = list1[0]

				hit = list1[2]
				if "NOG" in hit:
					pass
				else:
					hit = list1[3]

				if "." in hit:
					sub = hit.split(".")
					if "NOG" in hit:
						cog = sub[1]
					else:
						cog = sub[0]
				else:
					cog = hit


				bit_score = list1[5]
				score = float(bit_score)
						
				if score > 0:

					if score > bit_dict[query]:
						hit_dict[query] = cog
						bit_dict[query] = score
				else:
					print ids, hit, score

		bit_sorted = sorted(bit_dict.items(), key=operator.itemgetter(1), reverse=True)
		for item in bit_sorted:
			output.write(item[0] +"\t"+ str(hit_dict[item[0]]) +"\t"+ str(item[1]) +"\n")

		output.close()














