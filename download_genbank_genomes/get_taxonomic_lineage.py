import os
import sys
import re
from ete3 import NCBITaxa
ncbi = NCBITaxa()
from collections import defaultdict
import operator
import random
#ncbi.update_taxonomy_database()

refseq_main = open("/home/frankaylward/Desktop/marker_gene_benchmarking/assembly_summaries/Oct_20/final_set.txt", "r")

family_dict = defaultdict(list)
domain_dict = {}
phylum_dict = {}


for i in refseq_main.readlines():
	line = i.rstrip()
	if line.startswith("#"):
		pass
	else:
		list1 = line.split("\t")
		ident = list1[0]
		taxid = list1[9]
		goal = list1[13]
		ftp = list1[19]
		print taxid

		# lineage = list of taxids
		# name_list = list of names corresponding to the taxids

		lineage = ncbi.get_lineage(taxid)
		names = ncbi.get_taxid_translator(lineage)
		name_list = [names[taxid] for taxid in lineage]
		rank_dict = ncbi.get_rank(lineage)
		ranks = rank_dict.values()

		print lineage
		print names
		print rank_dict
		#print name_list
		print ""
		rank_str = " ".join(ranks)









