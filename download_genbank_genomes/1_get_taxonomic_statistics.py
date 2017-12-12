import os
import sys
import re
from ete3 import NCBITaxa
ncbi = NCBITaxa()
from collections import defaultdict
import operator
import random
#ncbi.update_taxonomy_database()

refseq_main = open("/home/frankaylward/Desktop/marker_gene_benchmarking/assembly_summaries/Oct_20/assembly_summary_genbank.txt", "r")
#refseq_hist = open("/home/frankaylward/Desktop/marker_gene_benchmarking/assembly_summaries/assembly_summary_refseq_historical.txt", "r")
output = open("genome_manifest.txt", "w")
subset = open("genome_manifest_partial.txt", "w")
exceptions = open("exceptions.txt", "w")
na_phyla = open("na_phyla.txt", "w")

family_tally = defaultdict(list)
tax = defaultdict(int)
c = defaultdict(int)
family_dict = defaultdict(list)
domain_dict = {}
phylum_dict = {}
manifest_dict = {}
tally = 0

for i in refseq_main.readlines():
	line = i.rstrip()
	if line.startswith("#"):
		pass
	else:
		list1 = line.split("\t")
		ident = list1[0]
		taxid = list1[5]
		goal = list1[13]
		ftp = list1[19]
		#tally +=1
		if ftp == "na":
			pass
		else:
			status = re.sub(" ", "_", list1[11])
			rest1 = "\t".join(list1[1:16])
			rest2 = "\t".join(list1[17:len(list1)])
			rest = rest1 +"\t"+ rest2
			#print rest2

			if goal == "Full":
				tally +=1
				try:
					lineage = ncbi.get_lineage(taxid)
					names = ncbi.get_taxid_translator(lineage)
					name_list = [names[taxid] for taxid in lineage]
					rank_dict = ncbi.get_rank(lineage)
					ranks = rank_dict.values()
					rank_str = " ".join(ranks)
					#print rank_str

					domain = str(name_list[2])
					tax[domain] +=1

					phylum = "NA"
					genus = "NA"
					superkingdom = "NA"

					if domain == "Bacteria" or domain == "Archaea" or domain == "Eukaryota":

						for m in rank_dict:
							#print rank_dict[m]
							if rank_dict[m] == "genus":
								genus = names[m]

							elif rank_dict[m] == "phylum":
								phylum = names[m]

							elif rank_dict[m] == "superkingdom":
								superkingdom = names[m]

						if phylum == "NA":
							na_phyla.write(ident +"\t"+ domain +"\t"+ superkingdom +"\t"+ phylum +"\t"+ genus +"\t"+ str(taxid) +"\t"+ list1[7] +"\n")
							if domain == "Eukaryota":
								pass
							else:
								output.write(ident +"\t"+ domain +"\t"+ superkingdom +"\t"+ phylum +"\t"+ genus +"\t"+ rest +"\n")
								subset.write(ident +"\t"+ domain +"\t"+ superkingdom +"\t"+ phylum +"\t"+ genus +"\t"+ str(taxid) +"\t"+ list1[7] +"\n")

						else:
							output.write(ident +"\t"+ domain +"\t"+ superkingdom +"\t"+ phylum +"\t"+ genus +"\t"+ rest +"\n")
							subset.write(ident +"\t"+ domain +"\t"+ superkingdom +"\t"+ phylum +"\t"+ genus +"\t"+ str(taxid) +"\t"+ list1[7] +"\n")

				except:
					exceptions.write(ident +"\t"+ rest +"\n")
					#print ident, taxid
print tally
# end















