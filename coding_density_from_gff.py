import sys, os
from collections import defaultdict
from Bio import SeqIO

fasta = open(sys.argv[1])
gff = open(sys.argv[2])

total_length = 0
coding_length = 0
status = 0

# first iterate through fasta file and get regions of N and the flanking regions
index_list = defaultdict(list)
for i in SeqIO.parse(fasta, "fasta"):
	seq = i.seq
	#print (seq)
	name = i.id
	for index,nucl in enumerate(list(seq)):
		if nucl == "N":
			index_list[i.id].append(index)

for i in index_list:
	new_values = []
	index = index_list[i]
	for j in index:
		prev = float(j) - 1
		next = float(j) + 1
		if prev in index and next in index:
			pass
			val +=1
		elif prev in index and not next in index:
			val +=1
			if val > 10:
				print (val)
				to_add_end = list(range(j, j+500))
				new_values = new_values + to_add_end
				to_add_start = list(range(startnum, startnum-500))
				new_values = new_values + to_add_start
				print(j, startnum)
			else:
				pass
			val = 0
		elif next in index and not prev in index:
			val = 0
			startnum = j
			#to_add = list(range(j-500, j))
			#new_values = new_values + to_add
			
	index_list[i] = set(index + new_values)
	#print (len(new_values), len(set(new_values)))
	#print(new_values)
	#print(index_list[i])

# now iterate through the gff file and get the coding density, excluding those regions of N
for i in gff.readlines():
	if status == 1:
		pass
	else:
		if i.startswith("##sequence-region"):
			line = i.rstrip()
			space = line.split(" ")
			name = space[1]
			length = float(space[3])
			effective_length = length - len(index_list[name])
			print (length, effective_length, len(index_list[name]))
			total_length += effective_length
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
				seq_name =  tabs[0]
				if start in index_list[seq_name] or end in index_list[seq_name]:
					pass
				else:
					diff = float(abs(end - start))
					coding_length += diff

if total_length < 200000:
	print("WARNING: Estimating coding density from only ", total_length)
print(coding_length, total_length, 100*(coding_length/total_length))