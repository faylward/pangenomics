import sys, os
from collections import defaultdict

input = open(sys.argv[1])

total_length = 0
coding_length = 0
status = 0

for i in input.readlines():
	if status == 1:
		pass
	else:
		if i.startswith("##sequence-region"):
			line = i.rstrip()
			space = line.split(" ")
			length = float(space[3])
			total_length += length
		elif i.startswith("##FASTA"):
			status = 1
		elif i.startswith("#"):
			pass
		else:
			line = i.rstrip()
			tabs = line.split("\t")
			if tabs[2] == "CDS" or tabs[2] == "rRNA" or tabs[2] == "tRNA" or tabs[2] == "tmRNA":
				#print(tabs)
				start = float(tabs[3])
				end = float(tabs[4])
				diff = float(abs(end - start))
				coding_length += diff
		
print(coding_length, total_length, 100*(coding_length/total_length))