import os
import sys
import subprocess
import re
import shlex
import operator
from collections import defaultdict
from Bio import SeqIO

sequences = open(sys.argv[1], "r")
cogs = open(sys.argv[2], "r")
outfolder = sys.argv[3]
trimvalue = sys.argv[4]

record_dict = SeqIO.to_dict(SeqIO.parse(sequences, "fasta"))

for i in cogs.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	
	example = tabs[0]
	underscore = example.split("_")
	cog = underscore[1]
	print cog
	filename = os.path.join(outfolder, cog+".faa")
	handle = open(filename, "w")

	records = [record_dict[j] for j in tabs]
	SeqIO.write(records, handle, "fasta")
	handle.close()

full_dict = {}
for i in os.listdir(outfolder):
	if i.endswith(".faa"):
		filename = os.path.join(outfolder, i)
		alignment = re.sub(".faa", ".aln", filename)
		trimmed = re.sub(".faa", ".trimal.aln", filename)

		cmd = "clustalo --force -i "+ filename +" -o "+ alignment 
		print cmd
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdout=open("output.txt", "w"), stderr=open("err.txt", "w"))

		cmd = "trimal -gt "+ trimvalue +" -in "+ alignment +" -out "+ trimmed
		print cmd
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdout=open("output.txt", "w"), stderr=open("err.txt", "w"))

		seq_dict = SeqIO.to_dict(SeqIO.parse(trimmed, "fasta"))
		full_dict.update(seq_dict)

for i in full_dict:
	print i, full_dict[i].id





















	



