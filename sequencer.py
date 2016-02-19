#!/usr/local/bin/python
'''
Created on Jan 8, 2016

@author: dulupinar
'''
import sequencingtools.tools as tools
import os.path
import cPickle as pickle
import logging as logger
import sys

#refReadFile = "./forcredit/hw1_W_2/ref_hw1_W_2_chr_1.txt"
#readFile = "./forcredit/hw1_W_2/reads_hw1_W_2_chr_1.txt"

refReadFile = "../HW1/practice_W_1/ref_practice_W_1_chr_1.txt"
readFile = "../HW1/practice_W_1/reads_practice_W_1_chr_1.txt"

refSeq = tools.ReferenceSequence(refReadFile)

if True:
    unmatchedReads = tools.readRead(readFile)
    kmerMap = tools.generateKmerMap(refSeq.refRead,50)
    i =0
    len_unmatchedReads = len(unmatchedReads)
    logger.warning("start sequencing")
    
    
    for read in unmatchedReads:
        #read is a tuple
        i+=1
        if i % (len_unmatchedReads/10) == 0:
            logger.warning("Process {0} has completed {1}% of read {2}".format(id,str(100*i/float(len_unmatchedReads)),1))
        
        found_read1 = refSeq.findMatch(read[0],kmerMap)
        if found_read1 == False:
            found_read1 = refSeq.findMatch(read[0],kmerMap,reverse = True)
                
        found_read2 = refSeq.findMatch(read[1],kmerMap)                       
        if found_read2 == False:
            found_read2 = refSeq.findMatch(read[1],kmerMap,reverse = True)
    
        
        #now we check for indels
        
        if (bool(found_read1) ^ bool(found_read2)):
            #only one was found
        
            if bool(found_read1) == True:
                #The pairedEnd had a match So lets check within a close distance
                start = found_read1
                check_read = read[1]
            else:
                start = found_read2
                check_read = read[0]
                
            logger.debug("checking for indel between {0} and {1}".format(str(start),str(start+400)))
            refSeq.findInDels(check_read,start + 50,start + 250) #need to check 50 after start
            logger.debug("done checking for indel")
   
    refSeq.printInfo()
