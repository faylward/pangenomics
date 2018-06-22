
import os, sys, re, subprocess, shlex
from collections import defaultdict
from Bio import SeqIO

input_folder = sys.argv[1]
outfolder = sys.argv[2]

for i in os.listdir(input_folder):
	if i.endswith(".fna") or i.endswith(".fa"):
		fasta = os.path.join(input_folder, i)
		core = re.sub(".fna", "", i)

		gff = os.path.join(outfolder, core+".output")
		prot = os.path.join(outfolder, core+".faa")
		nucl = os.path.join(outfolder, core+".genes.fna")

		cmd = "prodigal -i "+ fasta +" -f gff -o "+ gff +" -a "+ prot +" -d "+ nucl
		print i
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdout=open("prodigal.out", "w"), stderr=open("prodigal.err", "w"))

		gff_out = open(re.sub(".output", ".gff", gff), "w")
		gff_in = open(gff, "r")
		for j in gff_in.readlines():
			if j.startswith("#"):
				#gff_out.write(j)
				pass
			else:
				line = j.rstrip()
				tabs = line.split("\t")
				replicon = tabs[0]
				info = tabs[8]
				info_split = info.split(";")
				first = info_split[0]
				first_split = first.split("_")
				num = first_split[1]
				gene_id = replicon +"_"+ num

				last_info = ";".join(info_split[1:len(info_split)])
				new_info = "ID="+ gene_id +";" #+ last_info
				new_line = "\t".join(tabs[0:7]) +"\t.\t"+ new_info
				gff_out.write(new_line +"\n")


