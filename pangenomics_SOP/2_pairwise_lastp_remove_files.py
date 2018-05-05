
############################################################################################
##### This script computes pairwise homology between proteins in two genomes using LAST ####
##### It requires that LAST is installed and can be called from the working directory.  ####
##### Usage: pairwise_lastp.py <input_folder> <output_folder>                           ####
############################################################################################

import os, sys, re, subprocess, shlex
from collections import defaultdict
from Bio import SeqIO

input_folder = sys.argv[1]
output_folder = re.sub("_GENES", "_LASTP", input_folder)
length_dict = {}
proteins_in_genome = defaultdict(list)

cmd = "mkdir "+ output_folder
cmd2 = shlex.split(cmd)
subprocess.call(cmd2, stdin=open("stdout.txt", "w"), stdout=open("stderr.txt", "w"))

folders = os.listdir(input_folder)
print "Formatting LASTDBs..."
# First iterate through files and make lastdbs
for faa in folders:
	if faa.endswith(".faa"):
		prefix = re.sub(".faa", ".faa", faa)
		db = os.path.join(input_folder, prefix)
		filepath = os.path.join(input_folder, faa)

		# get length of each protein and put it in a dictionary. This is used for computing alingment lengths later. 
		for prot in SeqIO.parse(filepath, "fasta"):
			length_dict[prot.id] = len(prot.seq)
			proteins_in_genome[prefix].append(prot.id)

		# If the lastdb is already there, no need to compute it. 
		lastpath = os.path.join(input_folder, prefix+".suf")
		cmd = "lastdb -p "+ db +" "+ filepath
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdin=open("stdout.txt", "w"), stdout=open("stderr.txt", "w"))

###################################################################################
# Define function and dictionaries that will be used for parsing the last outputs #
###################################################################################
fullout = open("marinimicrobia.blast-graph", "w")
bit_dict1  = defaultdict(float)
eval_dict1 = defaultdict(float)
bit_dict2  = defaultdict(float)
eval_dict2 = defaultdict(float)
final_list = defaultdict(list)

def parse_lastout(genome1, genome2, lastout1, lastout2):

	forward = open(lastout1, "r")
	reverse = open(lastout2, "r")
	core = genome1 +"__"+ genome2

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

print "Running LASTAL searches..."
# Iterate through .faa files and make sure that the LAST index files have been computed. 
folders = os.listdir(input_folder)
# keep an ongoing tally of the comparisons
tally = 0
done = {}
for faa in folders:
	if faa.endswith(".faa"):
		#prefix = re.sub(".faa", "", faa)
		#db = os.path.join(input_folder, prefix)
		filepath = os.path.join(input_folder, faa)

		# Now iterate through all the files again and make the appropriate pairs. 
		for faa2 in folders:
			if faa2.endswith(".faa"):
				# No need to compare genomes to themselves (although that is sometimes done). 
				if faa == faa2:
					pass
				else:
					#prefix2 = re.sub(".faa", "", faa2)
					#db2 = os.path.join(input_folder, faa2)
					filepath2 = os.path.join(input_folder, faa2)

					if faa +"__"+ faa2 in done:
						pass
					else:
						# do the forward
						output1 = os.path.join(output_folder, faa +"__"+ faa2 +".lastout")
						cmd = "lastal -P 16 -m 50 "+ filepath2 +" "+ filepath +" -f BlastTab"
						cmd2 = shlex.split(cmd)
						subprocess.call(cmd2, stdin=open("stderr.txt", "w"), stdout=open(output1, "w"))

						# now do the reverse
						output2 = os.path.join(output_folder, faa2 +"__"+ faa +".lastout")
						cmd = "lastal -P 16 -m 50 "+ filepath +" "+ filepath2 +" -f BlastTab"
						#print cmd
						cmd2 = shlex.split(cmd)
						subprocess.call(cmd2, stdin=open("stderr.txt", "w"), stdout=open(output2, "w"))

						parse_lastout(faa, faa2, output1, output2)

						cmd = "rm "+ output1
						cmd2 = shlex.split(cmd)
						subprocess.call(cmd2, stdin=open("stderr.txt", "w"), stdout=open("stdout.txt", "w"))

						cmd = "rm "+ output2
						cmd2 = shlex.split(cmd)
						subprocess.call(cmd2, stdin=open("stderr.txt", "w"), stdout=open("stdout.txt", "w"))

						tally +=1
						print faa, faa2, tally
						done[faa +"__"+ faa2] = faa +"__"+ faa2
						done[faa2 +"__"+ faa] = faa2 +"__"+ faa

print "Compiling Results..."
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
		fullout.write("# "+first+"\t"+second+"\n")
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
				#print first, second
				proteins1 = proteins_in_genome[first]
				proteins2 = proteins_in_genome[second]

				if bit_dict1[pair] > 0 and bit_dict2[reverse] > 0:
					if query in proteins1:
						#print "worked"
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

# after this run proteinortho with the command: proteinortho5.pl -cpus=12 -step=3 -project=marinimicrobia














