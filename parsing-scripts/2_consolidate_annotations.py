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

cluster = open("/home/frankaylward/Desktop/marker_gene_benchmarking/Marinimicrobia/functional_annotations/longest_rep.txt", "r")

nog = open("/home/frankaylward/Desktop/marker_gene_benchmarking/Marinimicrobia/functional_annotations/hmmer_annotations/pangenome_vs_NOG.parsed", "r")
pfam = open("/home/frankaylward/Desktop/marker_gene_benchmarking/Marinimicrobia/functional_annotations/hmmer_annotations/pangenome_vs_pfam.parsed", "r")
tigr = open("/home/frankaylward/Desktop/marker_gene_benchmarking/Marinimicrobia/functional_annotations/hmmer_annotations/pangenome_vs_tigr.parsed", "r")

nog_annot = open("/home/frankaylward/Desktop/marker_gene_benchmarking/hmm/NOG/NOG.annotations.tsv", "r")
pfam_annot = open("/home/frankaylward/Desktop/marker_gene_benchmarking/hmm/Pfams/Pfam-A.clans.tsv", "r")
tigr_annot = open("/home/frankaylward/Desktop/marker_gene_benchmarking/hmm/TIGRfams/profile_annotations.tigr.txt", "r")

output = open("/home/frankaylward/Desktop/marker_gene_benchmarking/Marinimicrobia/functional_annotations/full_annotation.tsv", "w")
output.write("Protein_cluster\trepresentative\trep_length\tNOG\tNOG_score\tNOG_desc\tPfam\tPfam_score\tPfam_desc\tTIGR\tTIGR_score\tTIGR_desc\n")

################################################################
###### Loop through and parse the checkm HMM output ############
################################################################

# NOG
nog_desc = defaultdict(lambda:"NA")
for i in nog_annot.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	name = tabs[1]
	annot = tabs[5]
	nog_desc[name] = annot

# Pfam
pfam_desc = defaultdict(lambda:"NA")
for i in pfam_annot.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	name = tabs[0]
	annot = tabs[4]
	pfam_desc[name] = annot

# TIGR
tigr_desc = defaultdict(lambda:"NA")
for i in tigr_annot.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	name = tabs[0]
	annot = tabs[3]
	tigr_desc[name] = annot

################################################################
################################################################
################################################################

# get protein annotations
nog_dict = defaultdict(lambda:("NA", "NA"))
for i in nog.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	name = tabs[0]
	rest = (tabs[1], tabs[2])
	nog_dict[name] = rest

# get protein annotations
pfam_dict = defaultdict(lambda:("NA", "NA"))
for i in pfam.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	name = tabs[0]
	rest = (tabs[1], tabs[2])
	pfam_dict[name] = rest

# get protein annotations
tigr_dict = defaultdict(lambda:("NA", "NA"))
for i in tigr.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	name = tabs[0]
	rest = (tabs[1], tabs[2])
	tigr_dict[name] = rest

################################################################
################################################################
################################################################

for i in cluster.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	name = tabs[0]
	protein = tabs[1]
	length = tabs[2]

	nog = nog_dict[protein]
	pfam = pfam_dict[protein]
	tigr = tigr_dict[protein]

	#print pfam, protein

	nog_annot  =  nog_desc[nog[0]]
	pfam_annot = pfam_desc[pfam[0]]
	tigr_annot = tigr_desc[tigr[0]]

	print protein, nog

	output.write(name +"\t"+ protein +"\t"+ length +"\t"+ nog[0] +"\t"+ nog[1] +"\t"+ nog_annot +"\t"+ pfam[0] +"\t"+ pfam[1] +"\t"+ pfam_annot +"\t"+ tigr[0] +"\t"+ tigr[1] +"\t"+ tigr_annot +"\n")

	
























