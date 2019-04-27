#!/usr/bin/env python
import sys, os, re, shlex, subprocess, pandas, numpy, itertools, argparse, time
from collections import defaultdict
from Bio import SeqIO
import pandas

# predict proteins from genome FNA file
def predict_proteins(genome_file, project):

	seqdict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
	if len(seqdict) < 1:
		raise (genome_file+" does not appear to be in FASTA format!")

	file_name = os.path.basename(genome_file)
	file_base = os.path.splitext(file_name)[0]
	base_name = os.path.basename(file_base)
	protein_file = os.path.join(project, base_name+".faa")	
	
	cmd = "prodigal -i "+ genome_file +" -a "+ protein_file
	print(cmd)
	cmd2 = shlex.split(cmd)
	subprocess.call(cmd2, stdout=open("out.txt", "w"), stderr=open("err.txt", "w"))
	return protein_file
			
# run HMMER3
def run_hmmer(input_file, cpus, evalue):
	output_file = re.sub(".faa", ".hmmout", input_file)
	cmd = "hmmsearch -E "+ str(evalue) +" --cpu "+ cpus +" --tblout "+ output_file +" hmm/vog.small.hmm "+ input_file	
	print(cmd)
	cmd2 = shlex.split(cmd)
#	subprocess.call(cmd2, stdout=open("out.txt", "w"), stderr=open("err.txt", "w"))
	os.remove("out.txt")
	return output_file

# define function for parsing HMMER3 output
def parse_hmmout(hmmout):
	input = open(hmmout, "r")
	final_dict = defaultdict(int)
	hit_dict = defaultdict(lambda:"no_annot")
	bit_dict = defaultdict(float)

	for i in input.readlines():
		line = i.rstrip()
		if line.startswith("#"):
			pass
		else:
			newline = re.sub("\s+", "\t", line)
			tabs = newline.split("\t")
			protein = tabs[0]
			hit = tabs[2]
			eval = float(tabs[4])
			score = float(tabs[5])
			if score > bit_dict[protein]:
				bit_dict[protein] = score
				hit_dict[protein] = hit
			else:
				pass
	for i in hit_dict:
		vog = hit_dict[i]
		final_dict[vog] +=1
	return final_dict

df = pandas.DataFrame()
inputdir = sys.argv[1]
outputfile = sys.argv[2]
for i in os.listdir(inputdir):
	if i.endswith(".fna"):
		name = re.sub(".fna", "", i)
		inputfile = os.path.join(inputdir, i)
		print(inputfile)
		protein_file = predict_proteins(inputfile, inputdir)
		hmmout = run_hmmer(protein_file, "4", "1e-5")
		hit_dict = parse_hmmout(hmmout)
		
		s1 = pandas.DataFrame(pandas.Series(hit_dict, name = name))
		df = pandas.concat([df, s1], axis=1, sort=True)

df = df.fillna(0)
df.to_csv(outputfile, sep="\t")

