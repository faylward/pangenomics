
import os, sys, re, subprocess, shlex
from collections import defaultdict
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline

input1 = open(sys.argv[1], "r")
genes = open("all_genes.fna", "r")
proteins = open("all_proteins.faa", "r")
output_folder = "codon_alignments/"

cmd = "mkdir "+ output_folder
cmd2 = shlex.split(cmd)
subprocess.call(cmd2, stdin=open("stdout.txt", "w"), stdout=open("stderr.txt", "w"))

cmd = "mkdir "+ os.path.join(output_folder, "alignments")
cmd2 = shlex.split(cmd)
subprocess.call(cmd2, stdin=open("stdout.txt", "w"), stdout=open("stderr.txt", "w"))

gene_dict = SeqIO.to_dict(SeqIO.parse(genes, "fasta"))
protein_dict    = SeqIO.to_dict(SeqIO.parse(proteins, "fasta"))

clusters = ["cluster_00205", "cluster_00158", 'cluster_00013']

tally = 0
for i in input1.readlines():
	precord_list = []
	nrecord_list = []
	line = i.rstrip()
	tabs = line.split("\t")
	prots = tabs[9:len(tabs)]
	for j in protein_dict:
		if j in prots:
			precord = protein_dict[j]
			#seq = precord.seq
			#if "*" in seq:
			#	precord.seq = Seq("".join([x for x in seq if x != '*']))

			precord_list.append(precord)
			nrecord_list.append(gene_dict[j])

	name = tabs[0]
	if name in clusters:
		pass
		poutput = open(os.path.join(output_folder, "alignments", name+".faa"), "w")
		noutput = open(os.path.join(output_folder, "alignments", name+".fna"), "w")
		SeqIO.write(precord_list, poutput, "fasta")
		SeqIO.write(nrecord_list, noutput, "fasta")

poutput.close()
noutput.close()

for i in os.listdir(os.path.join(output_folder, "alignments")):
	if i.endswith(".faa"):
		faa = os.path.join(output_folder, "alignments", i)
		aa_aln = re.sub(".faa", ".faa.aln", faa)
		fna = re.sub(".faa", ".fna", faa)
		pal2nal = re.sub(".faa", ".pal2nal", faa)

		print faa
		#cmd = 'clustalo -t Protein --force --threads 6 -i '+ faa
		cmd = "clustalo -i "+ faa
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdout=open(aa_aln, 'w'), stderr=open('clustal.error.txt', 'w'),  env={'COLUMNS':'300'})

		cmd = 'pal2nal.pl '+ aa_aln +' '+ fna +' -output paml -nogap'
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdout=open(pal2nal, 'w'), stderr=open('error.txt', 'w'))

cmd = "cat "+ os.path.join("codon_alignments", "alignments", "*.pal2nal")
cmd2 = shlex.split(cmd)
subprocess.call(cmd, stdout=open("codon_alignments/all_paml.fna", "w"), stderr=open("stderr.txt", "w"), shell=True)


# after this go into the codon_alignments folder and type "codeml"



