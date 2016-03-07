'''
Created on Jan 8, 2016

@author: dulupinar
'''
import sequencingtools.tools as tools
import sys
import numpy as np
import os.path
import cPickle as pickle
import logging as logger
import operator
from operator import itemgetter
from _collections import defaultdict

def reject_outliers(data, m=2):
    return data[abs(data - np.mean(data)) < m * np.std(data)]


output = open("STR_FINAL_GRAD.txt",'w') 
refReadFile = "../HW2/practice_E_1/ref_practice_E_1_chr_1.txt"
readFile = "../HW2/practice_E_1/reads_practice_E_1_chr_1.txt"
#refReadFile = "../HW2/hw2grad/ref_hw2grad_M_1_chr_1.txt"

#create all the kmers
all_kmer = dict() #key is the dict the value is a class that contains all the information
tools.allKmerBetter(3, "", tools.alpha, all_kmer)
tools.allKmerBetter(4, "", tools.alpha, all_kmer)
tools.allKmerBetter(5, "", tools.alpha, all_kmer)
unmatchedReads = tools.readRead(readFile)
len_unmatchedReads = len(unmatchedReads)
i = 0


if True:
    #grad version
    refReadFile = "../HW2/hw2grad/ref_hw2grad_M_1_chr_1.txt"
    refSeq = tools.ReferenceSequence("donorReadin.txt")
    #refSeq.readIn("FINAL_GRAD.txt")
    donor = refSeq.generateDonor()
    refSeq.findSTRRegex() 
    refSeq.printInfo(output)
    sys.exit()
elif False:
    refSeq = tools.ReferenceSequence(refReadFile)  
    refSeq.readIn("str_test_grad_write_g_cov.1000.txt")
    donor = refSeq.generateDonor()
    refSeq.findSTRRegex() 
    refSeq.printInfo(output)
    sys.exit()
#####
refSeq = tools.ReferenceSequence(refReadFile)  
kmerMap = tools.generateKmerMap(refSeq.refRead,50)
unmatchedReads = tools.readRead(readFile)
i =0
len_unmatchedReads = len(unmatchedReads)


for read in unmatchedReads:
    #read is a tuple
    i+=1
    if i % (len_unmatchedReads/10) == 0:
        logger.warning("Process {0} has completed {1}% of read {2}".format(id,str(100*i/float(len_unmatchedReads)),1))
    
    found_read1 = refSeq.findMatch(read[0],kmerMap)
    read1_reverse = False
    if found_read1 == False:
        found_read1 = refSeq.findMatch(read[0],kmerMap,reverse = True)
        read1_reverse = True
           
    read2_reverse = False
    found_read2 = refSeq.findMatch(read[1],kmerMap)                       
    if found_read2 == False:
        found_read2 = refSeq.findMatch(read[1],kmerMap,reverse = True)
        read2_reverse = True

    
    #now check indel
    if (bool(found_read1) ^ bool(found_read2)):
        paired_reads = []
        if found_read1 == False:
            if read2_reverse == False:
                check_read = read[0][::-1]
            else:
                
                check_read = read[0]
            start = found_read2  
        if found_read2 == False:
            if read1_reverse == False:
                check_read = read[1][::-1]
            else:
                
                check_read = read[1]
            start = found_read1
        logger.debug("checking for indel between {0} and {1}".format(str(start),str(start+400)))
        refSeq.findInDels(check_read,start + 50,start + 250) #need to check 50 after start
        logger.debug("done checking for indel")

donor = refSeq.generateDonor()
refSeq.findSTRRegex() 
refSeq.printInfo(output)
