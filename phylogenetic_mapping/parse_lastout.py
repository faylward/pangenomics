import os, sys, re
from collections import defaultdict

input1 = open(sys.argv[1], "r")
bit_dict = defaultdict(float)
hit_dict = defaultdict(list)
percid_dict = defaultdict(float)
for i in input1.readlines():

	if i.startswith("#"):
		pass
	else:

		line = i.rstrip()
		tabs = line.split("\t")
		query = tabs[0]
		hit = tabs[1]
		bit = float(tabs[11])
		percid = tabs[2]

		if bit > 50:
			if bit >= bit_dict[query]:
				bit_dict[query] = float(bit)
				percid_dict[query] = percid

			if bit >= bit_dict[query]*0.95:
				hit_dict[query].append(hit)

for i in hit_dict:
	sys.stdout.write(i +"\t"+ str(percid_dict[i]) +"\t"+ str(bit_dict[i]) +"\t"+ ";".join(hit_dict[i]) +"\n")
