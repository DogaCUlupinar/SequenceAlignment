#!/usr/local/bin/python
'''
Created on Jan 8, 2016

@author: dulupinar
'''
import sequencingtools.tools as tools
import gc
import logging as logger
import glob
import sys
from __builtin__ import file

logger.basicConfig(level=logger.WARNING,format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

#refReadFile = "../HW2/practice_E_1/ref_practice_E_1_chr_1.txt"
#readFile = "../HW2/practice_E_1/reads_practice_E_1_chr_1.txt"

refReadFile = "../HW2/hw2grad/ref_hw2grad_M_1_chr_1.txt"
readPrefix = "../HW2/hw2grad/reads/splitreads*"
output = open("practice.txt",'w')  

readFiles = [file for file in glob.glob(readPrefix)]

refSeq = tools.ReferenceSequence(refReadFile)
logger.warning("reading in reference file: " + refReadFile)
kmerMap = tools.generateKmerMap(refSeq.refRead,50)
read_count = 0
for readFile in readFiles:
    logger.warning("reading in read file: " + readFile)
    unmatchedReads = tools.readRead(readFile)    
    len_unmatchedReads = len(unmatchedReads)
   
    logger.warning("have read in {0}% of files".format(str(readFiles.index(readFile))))
    i = 0
    logger.warning("starting sequencing")
    for read in unmatchedReads:
        #read is a tuple
        i+=1
        if i % (len_unmatchedReads/10) == 0:
            logger.warning("Process {0} has completed {1}% of read {2}".format(id,str(100*i/float(len_unmatchedReads)),read_count))
        
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
    read_count+=1
refSeq.printInfo(filestream=output)
