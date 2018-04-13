import os
import sys
import subprocess
import re
import shlex
from Bio import SeqIO

input_folder = "/home/frankaylward/Desktop/marker_gene_benchmarking/Marinimicrobia/all_genomes/genome_subsets_for_phylogeny_and_barnap/phylogeny/proteins"
output = open("/home/frankaylward/Desktop/marker_gene_benchmarking/Marinimicrobia/functional_annotations/protein_lengths.txt", "w")

for i in os.listdir(input_folder):
	if i.endswith(".faa"):
		finput = os.path.join(input_folder, i)
		for record in SeqIO.parse(finput, "fasta"):
			length = len(record.seq)
			output.write(record.id +"\t"+ str(length) +"\n")
