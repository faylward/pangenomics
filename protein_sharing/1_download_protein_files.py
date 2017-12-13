import os
import sys
import re
from collections import defaultdict
import operator
import subprocess
import shlex
#ncbi.update_taxonomy_database()

proteinortho = open("test_genome.proteinortho", "r")

for i in proteinortho.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	entries = tabs[3:len(tabs)]
	
