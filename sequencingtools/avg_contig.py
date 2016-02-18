#!/usr/local/bin/python
import re
filename = "output.txt"
pattern = re.compile("contig.*\s+([0-9]{1,3})$")
contig = []
with open(filename) as f:
	for line in f:
		matched = pattern.match(line)
		if matched:
			contig.append(int(matched.group(1)))
			print matched.group(1)
contig = sorted(contig)
print "The average is: ", str(sum(contig)/len(contig))
print "The mode is: ", str(contig[len(contig)/2])
print "The max is : ", str(max(contig))
