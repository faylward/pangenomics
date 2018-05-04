
import os, sys, re, subprocess, shlex
from collections import defaultdict
from Bio import SeqIO

input_folder = sys.argv[1]
outfolder = input_folder.rstrip("/")
outfolder = outfolder +"_GENES" 

cmd = "mkdir "+ outfolder
cmd2 = shlex.split(cmd)
subprocess.call(cmd2, stdin=open("stdout.txt", "w"), stdout=open("stderr.txt", "w"))

all_prot = []
all_nucl = []

for i in os.listdir(input_folder):
	if i.endswith(".fna") or i.endswith(".fa"):
		fasta = os.path.join(input_folder, i)
		core = re.sub(".fna", "", i)

		gff = os.path.join(outfolder, core+".output")
		prot = os.path.join(outfolder, core+".faa")
		nucl = os.path.join(outfolder, core+".genes.fna")

		# run Prodigal and predict the genes and their translations
		cmd = "prodigal -i "+ fasta +" -f gff -a "+ prot +" -d "+ nucl
		print i
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdout=open("prodigal.out", "w"), stderr=open("prodigal.err", "w"))

		# now iterate through the protein file and re-name the proteins so that the genome name can be easily retrieved.
		protein_names = open(re.sub(".faa", ".names.faa", prot), "w")
		gff = open(re.sub(".faa", ".gff", prot), "w")
		prot_in = open(prot, "r")
		prot_list = []
		for record in SeqIO.parse(prot, "fasta"):
			name = record.id
			newname = core +".."+ name
			record.id = newname
			prot_list.append(record)
			all_prot.append(record)

			desc = record.description
			space = desc.split(" ")

			if int(space[6]) == 1:
				sign = "+"
			else:
				sign = "-"

			gff.write(newname +"\tProdigal\tCDS\t"+ space[2] +"\t"+ space[4] +"\t0\t"+ sign +"\t.\tID="+ newname +";\n")

		# now write the new proteins to the file
		SeqIO.write(prot_list, protein_names, "fasta")

		# now iterate through the genes file and re-name the genes so that they are identical to those of the proteins
		gene_names = open(re.sub(".genes.fna", ".names.fna", nucl), "w")
		gene_in = open(nucl, "r")
		gene_list = []
		for record in SeqIO.parse(nucl, "fasta"):
			name = record.id
			newname = core +".."+ name
			record.id = newname
			gene_list.append(record)
			all_nucl.append(record)

		SeqIO.write(gene_list, gene_names, "fasta")

		# now delete the original protein and gene files to save space and avoid confusion in the future
		cmd = "rm "+ prot
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdout=open("stdout.txt", "w"), stderr=open("stderr.txt", "w"))

		cmd = "rm "+ nucl
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdout=open("stdout.txt", "w"), stderr=open("stderr.txt", "w"))

prot_output = open("all_proteins.faa", "w")
nucl_output = open("all_genes.fna", "w")

SeqIO.write(all_prot, prot_output, "fasta")
SeqIO.write(all_nucl, nucl_output, "fasta")
















