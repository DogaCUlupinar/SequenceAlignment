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
a = tools.generateKmerMap(refSeq.refRead,50)

for readFile in readFiles:
    logger.warning("reading in read file: " + readFile)
    unmatchedReads = tools.readRead(readFile)    
    len_unmatchedReads = len(unmatchedReads)
   
    logger.warning("have read in {0}% of files".format(str(readFiles.index(readFile))))
    i = 0
    logger.warning("starting sequencing")
    for read in unmatchedReads:
        i+=1
        if i % (len_unmatchedReads/25) == 0:
            logger.warning("we have completed " + str(100*i/float(len_unmatchedReads)) + "% of this read")
            
        if refSeq.findMatch(read,a) == False:
            if refSeq.findMatch(read,a,reverse = True) == False:
                if read.pairedRead.matchRangeInRef is not None:
    
                    #The pairedEnd had a match So lets check within a close distance
                    #start = read.pairedRead.matchRangeInRef[1]+50
                    #logger.debug("checking for indel between {0} and {1}".format(str(start),str(start+400)))
                    #refSeq.findInDels(read,start,start + 400)
                    logger.debug("done checking for indel")
                else:
                    #the pairedEnd did not find a match so lets check over the whole distance
                    #refSeq.findInDels(read,0,len(refSeq.refRead))
                    pass
         
refSeq.printInfo(filestream=output)
