import argparse
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

# This script requires that HMMER3 and Prodigal are install in the PATH. 

# create and wipe the LOGFILE before beginning
o = open("LOGFILE.txt", "w")
o.close()

#########################################
### define prodigal launcher function ###
#########################################
def predict_genes(folder):
	for filenames in os.listdir(folder):
		if filenames.endswith(".fna"):
			#print "Predicting genes..."
			input_file = os.path.join(folder, filenames)
			protein_file = re.sub(".fna", ".faa", input_file)
			gff_file = re.sub(".fna", ".gff", input_file)

			cmd = "prodigal -i "+ input_file +" -f gff -o "+ gff_file +" -a "+ protein_file 
			#print cmd
			cmd2 = shlex.split(cmd)
			subprocess.call(cmd2, stdout=open("LOGFILE.txt", 'w'), stderr=open("LOGFILE.txt", 'a'))
# end

#######################################
### define HMMER3 launcher function ###
#######################################
def run_hmmer(folder, db_path):
	for filenames in os.listdir(folder):
		if filenames.endswith(".faa"):
			input_file = os.path.join(folder, filenames)
			hmm_out = re.sub(".faa", ".hmm_out", input_file)
			dom_out = re.sub(".faa", ".hmm_dom_out", input_file)

			cmd = cmd = "hmmsearch --cpu 16 --tblout " + hmm_out + " --domtblout "+ dom_out + " " + db_path + " " + input_file
			cmd2 = shlex.split(cmd)
			subprocess.call(cmd2, stdout=open("LOGFILE.txt", 'w'), stderr=open("LOGFILE.txt", 'a'))
# end

######################################
######## get protein dictionary ######
######################################
def parse_faa(folder):
	seq_dict = {}
	for filenames in os.listdir(folder):
		if filenames.endswith(".faa"):
			input_file = os.path.join(folder, filenames)
			for j in SeqIO.parse(input_file, "fasta"):
				seq_dict[j.id] = j
	return seq_dict
# end

####################################
#### define HMM output parser ######
####################################
def hmm_parser(folder, suffix, output):
	cog_dict = defaultdict(list)
	score_list = {}
	prot_list = []
	combined_output = open(output, "w")
	#combined_output.write("protein\tacc\thit\tscore\tlength\tcategory\tspecies\tdomain\tphylum\torder\tp2\tp3\n")
	hits = []
	bit_dict = {}
	df = pandas.DataFrame()

	for filenames in os.listdir(folder):
		if filenames.endswith(suffix):

			acc = re.sub(suffix, "", filenames)
			f = open(folder+"/"+filenames, 'r')
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

					if score > 0:
						if score > bit_dict[ids]:
							hit_dict[ids] = hit
							bit_dict[ids] = score

			bit_sorted = sorted(bit_dict.items(), key=operator.itemgetter(1), reverse=True)
			output_list = []
			for item in bit_sorted:
				entry = item[0]
				output_list.append(entry +"\t"+ str(hit_dict[entry]) +"\t"+ str(bit_dict[entry]))

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
				#acc2 = taxon_dict2[acc]

				if nr in done:
					combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ score +"\tNH\t"+ "\n")
					#pass
				else:
					combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ score +"\tBH\t"+ "\n")
					done.append(nr)

			s1 = pandas.DataFrame(pandas.Series(hit_profile, name = acc))
			df = pandas.concat([df, s1], axis=1)
	return df, cog_dict
# end

cogs = ["COG0186", "COG0200", "COG0099", "COG0048", "COG0097", "COG0090", "COG0081", "COG0088", "COG0495", "COG0552", "COG0197", "COG0256", "COG0091", "COG0096", "COG0098", "COG0049", "COG0201", "COG0124", "COG0018", "COG0016", "COG0100", "COG0087", "COG0172", "COG0080", "COG0525", "COG0522", "COG0012", "COG0103", "COG0052", "COG0092", "COG0202", "COG0541", "COG0533", "COG0184", "COG0085", "COG0215", "COG0102", "COG0185", "COG0093", "COG0094"]



########################################################################
##### use argparse to run through the command line options given #######
########################################################################
args_parser = argparse.ArgumentParser(description="Script for identifying phylogenetic marker genes (PMGs) in sequencing data", epilog="Virginia Tech Department of Biological Sciences")
args_parser.add_argument('-i', '--input_folder', required=True, help='Input folder of FNA sequence files')
args_parser.add_argument('-o', '--output_folder', required=True, help='output folder.')
args_parser.add_argument('-d', '--db', required=True, help='Database for matching. Accepted sets are "embl", "checkm_bact", or "checkm_arch"')
args_parser = args_parser.parse_args()

# TODO: for other marker sets it will be necessary to add in a separate switch so that HMMER3 can run the --tc_cutoff option for checkm, rather than  using a separate file with cutoffs for the embl set. 

# set up object names for input/output/database folders
output_dir = args_parser.output_folder
input_dir = args_parser.input_folder
db = args_parser.db
db_path = os.path.join("hmm/embl/", db+".hmm")


###########################################
#############3## run functions ############
###########################################
print "Predicting genes..."
predict_genes(input_dir)

print "Running HMMER3..."
run_hmmer(input_dir, db_path)
speci_df, cog_dict = hmm_parser(input_dir, ".hmm_out", "all_hmm_out.txt")

# this is old code to output an annotation table- it may be useful in some cases but isn't part of the core code. 
#name2 = 'hmm_profile.speci.txt'
#speci_df.fillna(0, inplace=True, axis=1)
#df2 = speci_df.transpose()
#df2.to_csv(name2, sep='\t')

# get a dictionary that links protein IDs to their SeqIO records
seq_dict = parse_faa(input_dir)


###########################################################################
######### Final compilation of files for input to ete3 ####################
###########################################################################
hmm_out = open("all_hmm_out.txt", "r")
acc_dict = {}
cog_dict = defaultdict(list)
for j in hmm_out.readlines():
	line = j.rstrip()
	tabs = line.split("\t")
	hit_type = tabs[4]
	if hit_type == "BH":
		protein = tabs[0]
		acc = tabs[1]
		cog = tabs[2]
		acc_dict[protein] = acc
		#print protein, acc, cog
		cog_dict[cog].append(protein)

cog_file = hmm_out = open("cog_file.txt", "w")
out_proteins = open("proteins_for_phylogeny.faa", "w")
output_seqs = []
tally = []
tally2 = []
for n in cogs:
	cog = n
	protein_list = cog_dict[n]
	for i in protein_list:
		acc = acc_dict[i]
		acc2 = re.sub("_", ".", acc)
		new_name = acc2 +"_"+ cog
		cog_file.write(new_name +"\t")

		record = seq_dict[i]
		record.id = new_name
		record.description = ""
		output_seqs.append(record)

	cog_file.write("\n")

SeqIO.write(output_seqs, out_proteins, "fasta")



# and then run through ete3






