import os, sys, re
from Bio import Phylo
from collections import defaultdict

tree = Phylo.read('tree.nwk', 'newick')
#for i in tree.find_clades():
#	print i.support


input1 = open("lastout.parsed", "r")
node = defaultdict(float)
for i in input1.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	hits = tabs[3]
	hit_list = hits.split(";")
	mrca = tree.common_ancestor(hit_list)
	node[mrca] +=1
	
for i in node:
	print str(i) +"\t3\t"+ str(node[i]) +"\t#0000ff\t1\t0"
