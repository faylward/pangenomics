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
output = open("genome_manifest_2.txt", "w")
subset = open("genome_manifest_partial_2.txt", "w")
exceptions = open("exceptions.txt", "w")
na_phyla = open("na_phyla.txt", "w")

output.write("bioproject	biosample	wgs_master	refseq_category	taxid	species_taxid	organism_name	infraspecific_name	isolate	version_status	assembly_level	release_type	genome_rep	seq_rel_date	asm_name	gbrs_paired_asm	paired_asm_comp	ftp_path	excluded_from_refseq	relation_to_type_material	placement	cellular	domain	phylum	p1	p2	p3	p4	p5	p6	p7	p8	p9	p10	p11\n")

genus_tally = defaultdict(int)
len_list = []

for i in refseq_main.readlines():
	line = i.rstrip("\n")
	if line.startswith("#"):
		pass
	else:
		list1 = line.split("\t")
		ident = list1[0]
		taxid = list1[5]
		species_name = list1[7]
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
				try:
					lineage = ncbi.get_lineage(taxid)
					names = ncbi.get_taxid_translator(lineage)
					name_list = [names[taxid] for taxid in lineage]
					rank_dict = ncbi.get_rank(lineage)
					ranks = rank_dict.values()
					rank_str = " ".join(ranks)

					# lineage: list of taxids, each for a different taxonomic level
					# names: dictionary linking taxids to the names of the taxonomic levels 

					#print(lineage)
					#print(type(names))
					#print("  ".join(name_list))

					domain = str(name_list[2])

					if domain == "Bacteria" or domain == "Archaea": #or domain == "Eukaryota":

						slash = ftp.split("/")
						name2 = slash[9]

						#print("  ".join(name_list))
						for m in rank_dict:
							#print rank_dict[m]
							if rank_dict[m] == "species":
								#print rank_dict[m], m
								genus_tally[m] +=1
								if genus_tally[m] > 20:
									pass
								else:
									if len(name_list[0:7]) == 5:
										items = list(name_list[0:7])
										items.append("no_tax")
										items.append("no_tax")
										#print items
										name_item = "\t".join(items)
										#print name_item

										output.write(rest +"\t"+ species_name +"\t"+ name_item +"\n")
										subset.write(name2 +"\t"+ ident +"\t"+ species_name +"\t"+ name_item +"\n")

									elif len(name_list[0:7]) == 6:
										items = name_list[0:7]
										items.append("no_tax")
										name_item = "\t".join(items)
										#print items
										#print name_item

										output.write(rest +"\t"+ species_name +"\t"+ name_item +"\n")
										subset.write(name2 +"\t"+ ident +"\t"+ species_name +"\t"+  name_item +"\n")

									else:
										output.write(rest +"\t"+ species_name +"\t" +"\t".join(name_list[0:7]) +"\n")
										subset.write(name2 +"\t"+ ident +"\t"+ species_name +"\t" +"\t".join(name_list[0:7]) +"\n")
										#len_list.append(len(name_list[0:7]))


				except:
					exceptions.write(ident +"\t"+ rest +"\n")
					#print ident, taxid

# end















