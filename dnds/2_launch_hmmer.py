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

input_folder = sys.argv[1]
output_dir = sys.argv[2]
ncbi_file = open(sys.argv[3], "r")
speci_db = "/home/frankaylward/Desktop/marker_gene_benchmarking/hmm/speci/all.hmm"
info = open("/home/frankaylward/Desktop/marker_gene_benchmarking/hmm/speci/MG_BitScoreCutoffs.allhits.txt", "r")

# get taxonomic info
taxon_dict = {}
for j in ncbi_file:
	line = j.rstrip("\n")
	tabs = line.split("\t")
	ftp = tabs[19]
	slash = ftp.split("/")
	asm = slash[9]
	
	species = tabs[7]
	strain = tabs[8]
	name = species +"__"+ strain

	taxon_dict[asm] = name

#get cutoffs
markers = []
cutoff = {}
for i in info.readlines():
	if i.startswith("#"):
		pass
	else:
		line = i.rstrip()
		tabs = line.split("\t")
		name = tabs[0]
		cutoff[name] = float(tabs[1])
		markers.append(name)
		print name, tabs[1]

clear = open("error_file.txt", "w")
clear.close()

# define hmm launcher function
def hmm_launcher(folder):
	for files in os.listdir(input_folder):
		if files.endswith(".faa"):
			path = os.path.join(input_folder, files)
			tblout = os.path.join(output_dir, files+".hmmout")
			print path
	
			# run against the speci set
			cmd = "hmmsearch --cpu 16 --tblout " + tblout +" "+ speci_db + " " + path
			cmd2 = shlex.split(cmd)
			subprocess.call(cmd2, stdout=open("hmm.out", 'w'), stderr=open("error_file.txt", 'a'))

# end
hmm_launcher(input_folder)

################################################################
###### Loop through and parse the checkm HMM output ############
################################################################

def hmm_parser(folder, suffix, output):
	
	score_list = {}
	prot_list = []
	combined_output = open(output, "w")
	combined_output.write("protein\tacc\thit\tscore\tcategory\tspecies\n")
	hits = []
	bit_dict = {}
	df = pandas.DataFrame()

	for filenames in os.listdir(folder):
		if filenames.endswith(suffix):

			acc = re.sub(suffix, "", filenames)

			f = open(folder+"/"+filenames, 'r')
			o = open(folder+"/"+filenames+".parsed", 'w')
			hit_dict = {}
			bit_dict = defaultdict(int)
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
					if "COG" in hit:
						pass
					else:
						hit = acc2

					bit_score = list1[5]
					score = float(bit_score)

					if score > cutoff[hit]:

						if score > bit_dict[ids]:
							hit_dict[ids] = hit
							bit_dict[ids] = score
					else:
						print ids, hit, score

			bit_sorted = sorted(bit_dict.items(), key=operator.itemgetter(1), reverse=True)
			output_list = []
			for item in bit_sorted:
				entry = item[0]
				output_list.append(entry +"\t"+ str(hit_dict[entry]) +"\t"+ str(bit_dict[entry]))
				#o.write(entry +"\t"+ str(hit_dict[entry]) +"\t"+ str(bit_dict[entry]) +"\t"+ str(protein_lengths[entry]) +"\n")
			#o.close()

			#parsed = open(folder+"/"+filenames+".parsed", 'r')
			hit_profile = defaultdict(int)
			done = []
			for line in output_list:
				line1 = line.rstrip()
				tabs = line1.split("\t")
				ids = tabs[0]
				hits.append(ids)
				cog = tabs[1]
				score = tabs[2]
				nr = acc +"_"+ cog

				if nr in done:
					#combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ score +"\t"+ str(protein_lengths[ids]) +"\tNH\t"+ taxon_dict[acc] +"\n")
					#o.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ score +"\t"+ str(protein_lengths[ids]) +"\tNH\t"+ taxon_dict[acc] +"\n")
					pass
				else:
					combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ score +"\tBH\t"+ taxon_dict[acc] +"\n")
					#o.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ score +"\t"+ str(protein_lengths[ids]) +"\tBH\t"+ taxon_dict[acc] +"\n")
					done.append(nr)
			o.close()


	return df

#parse speci outputs
speci_df = hmm_parser(output_dir, "_genomic.faa.hmmout", "all_hmm_out.speci.txt")



