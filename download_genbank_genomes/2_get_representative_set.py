import os
import sys
import re
from ete3 import NCBITaxa
ncbi = NCBITaxa()
from collections import defaultdict
import operator
import random
#ncbi.update_taxonomy_database()

manifest = open("genome_manifest.txt", "r")
output = open("assembly_summary_subset.txt", "w")
final_set = open("final_set.txt", "w")

phylum_dict = defaultdict(list)
genus_dict = defaultdict(int)
manifest_dict = {}

for i in manifest.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	acc = tabs[0]
	rest = "\t".join(tabs[1:len(tabs)])
	manifest_dict[acc] = rest
	domain = tabs[1]
	phylum = tabs[3]
	genus = tabs[4]
	if domain == "Archaea" or domain == "Bacteria":
		if phylum == "NA":
			phylum_dict["unclassified_phylum"].append(acc)
		elif genus == "NA":
			phylum_dict[phylum].append(acc)
		else:
			if genus_dict[genus] > 20:
				pass
			else:
				genus_dict[genus] +=1
				phylum_dict[phylum].append(acc)
	else:
		phylum_dict["Eukaryotes"].append(acc)

for i in phylum_dict:
	length = len(phylum_dict[i])
	acc_list = ";".join(phylum_dict[i])
	output.write( i +"\t"+ str(length) +"\t"+ acc_list  +"\n" )

	if i == "unclassified_phylum":
		for j in phylum_dict[i]:
			final_set.write(j +"\t"+ manifest_dict[j] +"\n")
	else:
		accs = len(phylum_dict[i])
		if accs > 50:
			#print phylum_dict[i]
			subset = random.sample(phylum_dict[i], 50)
			for j in subset:
				final_set.write(j +"\t"+ manifest_dict[j] +"\n")
		else:
			for j in phylum_dict[i]:
				final_set.write(j +"\t"+ manifest_dict[j] +"\n")

