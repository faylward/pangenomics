
############################################################################################
##### This script computes pairwise homology between proteins in two genomes using LAST ####
##### It requires that LAST is installed and can be called from the working directory.  ####
##### Usage: pairwise_lastp.py <input_folder> <output_folder>                           ####
############################################################################################

import os, sys, re, subprocess, shlex
from collections import defaultdict
from Bio import SeqIO
import numpy as np

input_folder = sys.argv[1]
output_folder = sys.argv[2]
length_dict = {}
proteins_in_genome = defaultdict(list)

# define lastout parsing function
def parse_lastout(lastout):
		bit_dict = defaultdict(float)
		hit_dict = defaultdict(list)
		perc_dict = defaultdict(float)
		done = {}
		
		handle = open(lastout, "r")
		for j in handle.readlines():
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

#				if query in done:
#					pass
#				else:
				done[query] = query
				if evalue < 1e-3:
					if bit >= bit_dict[query]:
						hit_dict[query].append(hit)
						bit_dict[query] = bit
						perc_dict[query] = perc
					else:
						pass
						
		return hit_dict, bit_dict, perc_dict


# Iterate through .faa files and make sure that the LAST index files have been computed. 
folders = os.listdir(input_folder)
print "Formatting LAST databases..."
for faa in folders:
	if faa.endswith(".faa"):
		prefix = re.sub(".faa", "", faa)
		db = os.path.join(input_folder, prefix)
		filepath = os.path.join(input_folder, faa)

		# get length of each protein and put it in a dictionary. This is used for computing alingment lengths later. 
		for prot in SeqIO.parse(filepath, "fasta"):
			length_dict[prot.id] = len(prot.seq)
			proteins_in_genome[prefix].append(prot.id)

		lastpath = os.path.join(input_folder, prefix+".suf")
		cmd = "lastdb -p "+ db +" "+ filepath
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdin=open("stdout.txt", "w"), stdout=open("stderr.txt", "w"))

done = {}
for faa in folders:
	if faa.endswith(".faa"):
		prefix = re.sub(".faa", "", faa)
		db = os.path.join(input_folder, prefix)
		filepath = os.path.join(input_folder, faa)
		
		for faa2 in folders:
		
			if faa2.endswith(".faa"):
				prefix2 = re.sub(".faa", "", faa2)
				db2 = os.path.join(input_folder, prefix2)
				filepath2 = os.path.join(input_folder, faa2)

				if prefix +"__"+ prefix2 in done:
					pass
				elif prefix == prefix2:
					pass
				else:
					done[prefix +"__"+ prefix2] = prefix +"__"+ prefix2
					done[prefix2 +"__"+ prefix] = prefix2 +"__"+ prefix

					# run first lastout
					output1 = os.path.join(output_folder, prefix +"__"+ prefix2 +".lastout")
					cmd = "lastal -P 8 -m 500 "+ db2 +" "+ filepath +" -f BlastTab"
					cmd2 = shlex.split(cmd)
				#	subprocess.call(cmd2, stdin=open("stderr.txt", "w"), stdout=open(output1, "w"))

					# now run reciprocal
					output2 = os.path.join(output_folder, prefix2 +"__"+ prefix +".lastout")
					cmd = "lastal -P 8 -m 500 "+ db +" "+ filepath2 +" -f BlastTab"
					cmd2 = shlex.split(cmd)
				#	subprocess.call(cmd2, stdin=open("stderr.txt", "w"), stdout=open(output2, "w"))

					hit_dict1, bit_dict1, perc_dict1 = parse_lastout(output1)
					hit_dict2, bit_dict2, perc_dict2 = parse_lastout(output2)

					aai_dict = defaultdict(list)

					for query1 in hit_dict1:
						hit_list1 = hit_dict1[query1]
						
						hit_list2 = []
						for hit1 in hit_list1:
							if query1 in hit_dict2[hit1]:
						
							#if hit1 in hit_dict2:
								#print query1, hit1
								#hit2 = hit_dict2[hit1]
							#	query2 = hit1
							#	hit2 = hit_dict2[query2]
								#print hit1, hit2, query1, query2
							#	if query1 == hit2 and query2 == hit1:
								meanval = np.mean([perc_dict1[query1], perc_dict2[hit1]])
								aai_dict[prefix +"__"+ prefix2].append(meanval)
								break
							#else:
							#	print query1, query2, hit1, hit2
								
					aai = np.mean(aai_dict[prefix +"__"+ prefix2])
					num_hits = float(len(aai_dict[prefix +"__"+ prefix2]))
					num_prots1 = float(len(proteins_in_genome[prefix]))
					num_prots2 = float(len(proteins_in_genome[prefix2]))
					af1 = 100*(num_hits / num_prots1)
					af2 = 100*(num_hits / num_prots2)
					
					print prefix, prefix2, aai, num_hits, num_prots1, num_prots2, af1, af2



