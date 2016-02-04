#!/usr/local/bin/python
from collections import *

def problem1(a,b):
	#return the sum of all odd numbers from a to b
	sum = 0 
	for i in range(a,b+1):
		if (i % 2 == 1):
			sum+=i
	return sum

def problem2(filename):
	#return all the even numbers from a file
	linenum = 1
	str_content = ""
	with open(filename,'r') as f:
		for line in f:
			if linenum % 2 == 0:
				str_content+=line
			linenum+=1
	return str_content.strip()

def problem3(string):
	#counting the number of time each word appears in string
	word_list = string.split()
	cnt = Counter()
	for word in word_list:
		cnt[word]+=1

	return_string = ""
	for key,value in cnt.iteritems():
		return_string+="{0} {1}\n".format(key,value)
	return return_string.strip()

def problem4(string):
	#counting the number of each base in a strand of dna (string)
	cnt = Counter()
	for letter in string:
		cnt[letter]+=1

	return_string = "{0} {1} {2} {3}".format(cnt["A"],cnt['C'],cnt['G'],cnt['T'])
	return return_string.strip()

def problem5(dna,motif):
	#finding motifs in dna
	len_motif = len(motif)
	index = 0
	dd = defaultdict(list)
	while index + len_motif <= len(dna):		
		dd[dna[index:index+len_motif]].append(index + 1)
		index+=1

	return " ".join(map(str, dd[motif]))



if __name__ == "__main__":
	dna = "ACGAATACCCCGAATACCGAATACCGAATACCGGCGAATACACGAATACGTCGGACGAATACATACGAATACAACACACGCGAATACTAGCGAATACCGAATACAATAGCGAATACACGAATACAGACCGAATACAGACGAATACAACCGGTCGAATACCGAATACATGACCGAATACTAATCGAATACGGCTCAAATTATACCGAATACCGCGAATACGCCGAATACGCCATCAGTCCCGAATACGTTCGAATACCCGAATACCTATCGAATACACGAATACTATGTTCGAATACATAGCTAACGAATACCTTCGAATACCCCATTTCCGAATACGCGAATACTGCGGGCGTTCCGAATACCCATGATATCGCGAATACCACTCGAATACCGAATACGCCGAATACCGAATACCGAATACACGAATACCGACTTCGACGAATACCGAATACACGAATACCGCACGAATACTTGTTGCGAATACGTCCCTCTTCGAATACGACGAATACCGAATACCCTCGGTACGAATACAAAATGTGACGAATACCTCGAATACGACGAATACCGAATACCGAATACTCGAATACGTCCGAATACCGAATACTAACTACGAATACCATCGCTTAGCACCGAATACACGAATACCGAATACCCGAATACCGAATACTGCGAATACGCAGGCTAGCCGAATACATCGAATACTTCGCGAATACCTTCGAATACCCACGAATACAGGCGCGAATACCGAATACCCTGGGCGAATACCGAATACTCCGAATACACGAATACCTGTCGAATACCCGAATACACGAAGTCGAATACGACGAATACATTCCGAATACCGAATACAAAGCAAAAGCCATCGAATAC"
	motif = "CGAATACCG"
	print problem5(dna,motif)