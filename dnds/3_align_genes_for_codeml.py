
import os, sys, re, subprocess, shlex
from collections import defaultdict
from Bio import SeqIO

hmmout = open(sys.argv[1], "r")
genes = sys.argv[2]
output_folder = "data_products/"

name2cog = {}
for i in hmmout.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	name = tabs[0]
	cog = tabs[2]
	if cog == "COG0086":
		pass
	else:
		name2cog[name] = cog

cog2proteins = defaultdict(list)
for files in os.listdir(genes):
	if files.endswith(".faa"):
		faa = os.path.join(genes, files)
		for record in SeqIO.parse(faa, "fasta"):
			if record.id in name2cog:
				cog = name2cog[record.id]
				cog2proteins[cog].append(record)
	
cog2genes = defaultdict(list)
for files in os.listdir(genes):
	if files.endswith(".genes.fna"):
		fna = os.path.join(genes, files)
		for record in SeqIO.parse(fna, "fasta"):
			if record.id in name2cog:
				cog = name2cog[record.id]
				cog2genes[cog].append(record)

for item in cog2proteins:
	output = open(os.path.join(output_folder, item+".faa"), "w")
	SeqIO.write(cog2proteins[item], output, "fasta")

for item in cog2genes:
	output = open(os.path.join(output_folder, item+".fna"), "w")
	SeqIO.write(cog2genes[item], output, "fasta")

for i in os.listdir(output_folder):
	if i.endswith(".faa"):
		faa = os.path.join(output_folder, i)
		aa_aln = re.sub(".faa", ".faa.aln", faa)
		fna = re.sub(".faa", ".fna", faa)
		pal2nal = re.sub(".faa", ".pal2nal", faa)

		cmd = 'clustalo -i '+ faa +' -o '+ aa_aln
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdout=open('output.txt', 'w'), stderr=open('error.txt', 'w'))

		cmd = 'pal2nal.pl '+ aa_aln +' '+ fna +' -output paml -nogap'
		cmd2 = shlex.split(cmd)
		print cmd
		print pal2nal
		subprocess.call(cmd2, stdout=open(pal2nal, 'w'), stderr=open('error.txt', 'w'))









