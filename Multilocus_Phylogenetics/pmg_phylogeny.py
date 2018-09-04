#!/home/frankaylward/anaconda_ete/bin/python
import argparse, os, sys, subprocess, re, shlex, pandas, operator
from collections import defaultdict
from Bio import SeqIO

# specify directory where the script is located
homedir = "/home/frankaylward/Desktop/marker_gene_benchmarking/pangenomics/Multilocus_Phylogenetics"

# This script requires that HMMER3 and Prodigal are installed in the PATH. 

# create and wipe the LOGFILE before beginning
o = open("LOGFILE.txt", "w")
o.close()

#########################################
### define prodigal launcher function ###
#########################################
suffix = ".fa"
def predict_genes(folder):
	for filenames in os.listdir(folder):
		if filenames.endswith(".fna") or filenames.endswith(".fa"):
			s = filenames.split(".")
			suffix = "."+s[-1]

			#print "Predicting genes..."
			input_file = os.path.join(folder, filenames)
			protein_file = re.sub(suffix, ".faa", input_file)
			gff_file = re.sub(suffix, ".gff", input_file)
			acc = re.sub(suffix, "", filenames)

			cmd = "prodigal -i "+ input_file +" -f gff -o "+ gff_file +" -a "+ protein_file 
			#print cmd
			cmd2 = shlex.split(cmd)
			subprocess.call(cmd2, stdout=open("LOGFILE.txt", 'w'), stderr=open("LOGFILE.txt", 'a'))


def make_nr(folder):
	for filenames in os.listdir(folder):
		if filenames.endswith(".faa"):

			protein_file = os.path.join(folder, filenames)
			acc1 = re.sub(".faa", "", filenames)
			acc = re.sub("_protein", "", acc1)

			# make protein names non-redundant
			output_seqs = []
			handle = open(protein_file, "r")
			for record in SeqIO.parse(handle, "fasta"):
				name = record.id
				if ".." in name:
					pass
				else:
					record.id = acc +".."+ name
					#print record.id, acc, filenames, name

				output_seqs.append(record)

			handle.close()
			output = open(protein_file, "w")
			SeqIO.write(output_seqs, output, "fasta")

# end

#######################################
### define HMMER3 launcher function ###
#######################################
def run_hmmer(folder, db_path, cutoff):
	for filenames in os.listdir(folder):
		if filenames.endswith(".faa"):
			input_file = os.path.join(folder, filenames)
			hmm_out = re.sub(".faa", ".hmm_out", input_file)
			dom_out = re.sub(".faa", ".hmm_dom_out", input_file)

			cmd = cmd = "hmmsearch --cpu 16 --tblout " + hmm_out + " "+ cutoff +" --domtblout "+ dom_out + " " + db_path + " " + input_file
			cmd2 = shlex.split(cmd)
			subprocess.call(cmd2, stdout=open("LOGFILE.txt", 'w'), stderr=open("LOGFILE.txt", 'a'))
# end

######################################
######## get protein dictionary ######
######################################
def parse_faa(folder, dictionary):
	seq_dict = {}
	for filenames in os.listdir(folder):
		if filenames.endswith(".faa"):
			input_file = os.path.join(folder, filenames)
			for j in SeqIO.parse(input_file, "fasta"):
				#print j.id
				if j.id in dictionary:
					#print j.id
					seq_dict[j.id] = j
	return seq_dict
# end

####################################
#### define HMM output parser ######
####################################
def hmm_parser(folder, suffix, output):
	full_hit_dict = {}
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
			#bit_dict = {}
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

					if score > bit_dict[ids]:
						hit_dict[ids] = hit
						bit_dict[ids] = score

			bit_sorted = sorted(bit_dict.items(), key=operator.itemgetter(1), reverse=True)
			output_list = []
			for item in bit_sorted:
				entry = item[0]
				if entry in hit_dict:
					#if entry in bit_dict and entry in hit_dict:
					output_list.append(entry +"\t"+ str(hit_dict[entry]) +"\t"+ str(bit_dict[entry]))
					full_hit_dict[entry] = entry
					#print entry
				#else:
				#	print entry

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
					combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ score +"\tNH\t"+ "\n")
					#pass
				else:
					combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ score +"\tBH\t"+ "\n")
					done.append(nr)

			#s1 = pandas.DataFrame(pandas.Series(hit_profile, name = acc))
			#df = pandas.concat([df, s1], axis=1)
	return full_hit_dict
# end

########################################################################
##### use argparse to run through the command line options given #######
########################################################################
args_parser = argparse.ArgumentParser(description="Script for identifying phylogenetic marker genes (PMGs) in sequencing data\nFrank O. Aylward, Assistant Professor, Virginia Tech Department of Biological Sciences <faylward at vt dot edu>", epilog="Virginia Tech Department of Biological Sciences")
args_parser.add_argument('-i', '--input_folder', required=True, help='Input folder of FNA sequence files')
args_parser.add_argument('-c', '--cogs', required=True, help='output file for the COGs file, used in ete3')
args_parser.add_argument('-o', '--output_faa', required=True, help='amino acid fasta output file of the PMGs identified')
args_parser.add_argument('-f', '--genes', type=bool, default=False, const=True, nargs='?', help='Implies genes have already been predicted, and .faa files are already in the input folder)')
args_parser.add_argument('-n', '--nr', type=bool, default=False, const=True, nargs='?', help='Implies genes are already non-redundant and no re-naming is necessary)')
args_parser.add_argument('-d', '--db', required=True, help='Database for matching. Accepted sets are "embl", "checkm_bact", "checkm_arch", "rpl1_arch", "rpl1_bact", or "rpl2"')
args_parser = args_parser.parse_args()

# set up object names for input/output/database folders
input_dir = args_parser.input_folder
cogs_file = args_parser.cogs
faa_out = args_parser.output_faa
db = args_parser.db
db_path = os.path.join(homedir, "hmm/", db+".hmm")

######## Check if genes have already been predicted.
if not (args_parser.genes):
	print "Predicting genes...\n"
	predict_genes(input_dir)
else:
	print "Genes already predicted.\n"

####### Check if predicted proteins should be made non-redundant
if not (args_parser.nr):
	print "Making protein entries non-redundant...\n"
	make_nr(input_dir)
else:
	print "Genes already non-redundant.\n"

###### Get a list of the HMM entries that should be used
cogs = []
list_file = open(os.path.join(homedir, "hmm/", db+".list"), "r")
for i in list_file.readlines():
	line = i.rstrip()
	cogs.append(line)

###########################################
################ run functions ############
###########################################
cutoff_dict = {}
if "checkm" in db:
	cutoff = "--cut_nc"
else:
	cutoff = " "
	cutoff_file = open(os.path.join(homedir, "hmm/embl_cutoffs.txt"), "r")
	for i in cutoff_file.readlines():
		line = i.rstrip()
		tabs = line.split("\t")
		name = tabs[0]
		value = tabs[1]
		cutoff_dict[name] = float(value)

print "Running HMMER3...\n"
run_hmmer(input_dir, db_path, cutoff)

print "Parsing HMMER3 outputs...\n"
all_hits = hmm_parser(input_dir, ".hmm_out", "all_hmm_out.txt")
#print all_hits

# this is old code to output an annotation table- it may be useful in some cases but isn't part of the core code. 
#name2 = 'hmm_profile.speci.txt'
#df.fillna(0, inplace=True, axis=1)
#df2 = speci_df.transpose()
#df2.to_csv(name2, sep='\t')

print "Getting protein dictionary...\n"
# get a dictionary that links protein IDs to their SeqIO records
seq_dict = parse_faa(input_dir, all_hits)

###########################################################################
######### Final compilation of files for input to ete3 ####################
###########################################################################
print "Final compilation...\n"
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
		score = tabs[3]
		if db == "embl":
			if score >= cutoff_dict[cog]:
				acc_dict[protein] = acc
				#print protein, acc, cog
				cog_dict[cog].append(protein)
		else:
			acc_dict[protein] = acc
			#print protein, acc, cog
			cog_dict[cog].append(protein)
#print cog_dict

cog_file = hmm_out = open(cogs_file, "w")
out_proteins = open(faa_out, "w")
output_seqs = []
tally = []
tally2 = []
already_done = {}
for n in cogs:
	protein_list = cog_dict[n]
	#print protein_list
	for i in protein_list:
		acc = acc_dict[i]
		acc2 = re.sub("_", ".", acc)
		new_name = acc2 +"_"+ n
		
		if new_name in already_done:
			print n, i, acc, acc2

		else:
			already_done[new_name] = new_name

			#if n == "TIGR01128" and acc == "Marinimicrobia_Bin_45-1":
			#	print "Original\t", n, i, acc, acc2

			cog_file.write(new_name +"\t")
			tally.append(new_name)

			record = seq_dict[i]
			record.description = record.id
			#record.description = ""
			record.id = new_name
			output_seqs.append(record)
			tally2.append(record.id)

	cog_file.write("\n")

SeqIO.write(output_seqs, out_proteins, "fasta")
print len(tally), len(tally2), len(set(tally)), len(set(tally2))


# and then run through ete3 using the command "ete3 build -a proteins_for_phylogeny.faa --cogs cog_file.txt -w standard_trimmed_fasttree -m sptree_fasttree_90 --clearall -C 8 -o test_tree"






