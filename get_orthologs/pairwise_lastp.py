
############################################################################################
##### This script computes pairwise homology between proteins in two genomes using LAST ####
##### It requires that LAST is installed and can be called from the working directory.  ####
##### Usage: pairwise_lastp.py <input_folder> <output_folder>                           ####
############################################################################################

import os, sys, re, subprocess, shlex
from collections import defaultdict
from Bio import SeqIO

input_folder = sys.argv[1]
output_folder = sys.argv[2]
length_dict = {}
proteins_in_genome = defaultdict(list)

# Iterate through .faa files and make sure that the LAST index files have been computed. 
folders = os.listdir(input_folder)
for faa in folders:
	if faa.endswith(".faa"):
		prefix = re.sub(".faa", "", faa)
		db = os.path.join(input_folder, prefix)
		filepath = os.path.join(input_folder, faa)

		# get length of each protein and put it in a dictionary. This is used for computing alingment lengths later. 
		for prot in SeqIO.parse(filepath, "fasta"):
			length_dict[prot.id] = len(prot.seq)
			proteins_in_genome[prefix].append(prot.id)

		# If the lastdb is already there, no need to compute it. 
		lastpath = os.path.join(input_folder, prefix+".suf")
		if os.path.isfile(lastpath):
			pass
		else:
			cmd = "lastdb -p "+ db +" "+ filepath
			print cmd
			cmd2 = shlex.split(cmd)
			subprocess.call(cmd2, stdin=open("stdout.txt", "w"), stdout=open("stderr.txt", "w"))

		# Now iterate through all the files again and make the appropriate pairs. 
		for faa2 in folders:
			if faa2.endswith(".faa"):
				# No need to compare genomes to themselves (although that is sometimes done). 
				if faa == faa2:
					pass
				else:
					prefix2 = re.sub(".faa", "", faa2)
					db2 = os.path.join(input_folder, prefix2)
					filepath2 = os.path.join(input_folder, faa2)
					lastpath2 = os.path.join(input_folder, prefix2+".suf")
					if os.path.isfile(lastpath):
						pass
					else:
						cmd = "lastdb -p "+ db2 +" "+ filepath2
						cmd2 = shlex.split(cmd)
						subprocess.call(cmd2, stdin=open("stdout.txt", "w"), stdout=open("stderr.txt", "w"))

					# Now that lastdbs should be available for both the query and the reference, we can run the last command. 
					output = os.path.join(output_folder, prefix +"__"+ prefix2 +".lastout")
					cmd = "lastal -P 8 -m 500 "+ db2 +" "+ filepath +" -f BlastTab"
					print cmd
					cmd2 = shlex.split(cmd)
					subprocess.call(cmd2, stdin=open("stderr.txt", "w"), stdout=open(output, "w"))

# Now that the LAST comparisons have been computed we need to parse through the outputs and report the best hits to a consolidated file
fullout = open("marinimicrobia.blast-graph", "w")
bit_dict1  = defaultdict(float)
eval_dict1 = defaultdict(float)
bit_dict2  = defaultdict(float)
eval_dict2 = defaultdict(float)
final_list = defaultdict(list)

for i in os.listdir(output_folder):
	if i.endswith(".lastout"):

		forward = open(os.path.join(output_folder, i))
		core = re.sub(".lastout", "", i)
		pieces = core.split("__")
		print core
		first = pieces[0]
		second = pieces[1]
		reverse_core = second +"__"+ first +".lastout"
		reverse = open(os.path.join(output_folder, reverse_core))

		# iterate through the forward lastout
		bit_score = defaultdict(float)
		for j in forward.readlines():
			if j.startswith("#"):
				pass
			else:
				line = j.rstrip()
				tabs = line.split("\t")
				query = tabs[0]
				hit = tabs[1]
				perc = float(tabs[2])
				evalue = float(tabs[10])
				bit = float(tabs[11])
				aln = float(tabs[3])
				cov1 = aln / length_dict[hit]
				cov2 = aln / length_dict[query]
				queryhit  = query +"__"+ hit

				if evalue < 1e-5 and float(cov1) > 0.5 and bit > bit_score[query]*0.95:
					bit_dict1[queryhit] = bit
					eval_dict1[queryhit] = evalue

					if bit > bit_score[query]:
						bit_score[query] = bit
						final_list[core].append(queryhit)

		# iterate through the reverse lastout
		bit_score = defaultdict(float)
		for j in reverse.readlines():
			if j.startswith("#"):
				pass
			else:
				line = j.rstrip()
				tabs = line.split("\t")
				query = tabs[0]
				hit = tabs[1]
				perc = float(tabs[2])
				evalue = float(tabs[10])
				bit = float(tabs[11])
				aln = float(tabs[3])
				cov1 = aln / length_dict[hit]
				cov2 = aln / length_dict[query]
				queryhit = query +"__"+ hit

				if evalue < 1e-5 and float(cov1) > 0.5 and bit > bit_score[query]*0.95:
					bit_dict2[queryhit] = bit
					eval_dict2[queryhit] = evalue

					if bit > bit_score[query]:
						bit_score[query] = bit
						final_list[core].append(queryhit)

# Now that we have all of the hits in the appropriate dictionaries, we can consolidate and produce a final blast-graph output for proteinortho to use. 
done = {}
for lists in final_list:
	pieces = lists.split("__")
	first = pieces[0]
	second = pieces[1]
	print first, second
	if first +"__"+ second in done or second +"__"+ first in done:
	#if 1==0:	
		pass
	else:
		fullout.write("# "+first+".faa\t"+second+".faa\n")
		done[first +"__"+ second] = first +"__"+ second

		allhits = final_list[lists]
		for pair in allhits:
			pie = pair.split("__")
			query = pie[0]
			hit = pie[1]
			reverse = hit +"__"+ query

			if pair in done or reverse in done:
				pass
			else:
				done[pair] = pair
				done[reverse] = reverse
				proteins1 = proteins_in_genome[first]
				proteins2 = proteins_in_genome[second]

				if bit_dict1[pair] > 0 and bit_dict2[reverse] > 0:
					if query in proteins1:
						bit1 = bit_dict1[pair]
						bit2 = bit_dict2[reverse]
						eval1 = eval_dict1[pair]
						eval2 = eval_dict2[reverse]
						fullout.write(query +"\t"+ hit +"\t"+ str(eval1) +"\t"+ str(bit1) +"\t"+ str(eval2) +"\t"+ str(bit2) +"\n")

					elif query in proteins2:
						bit1 = bit_dict1[pair]
						bit2 = bit_dict2[reverse]
						eval1 = eval_dict1[pair]
						eval2 = eval_dict2[reverse]
						fullout.write(hit +"\t"+ query +"\t"+ str(eval2) +"\t"+ str(bit2) +"\t"+ str(eval1) +"\t"+ str(bit1) +"\n")
	
					else:
						print query
















