'''
Created on Jan 8, 2016

@author: dulupinar
'''
import sequencingtools.tools as tools
import os.path
import cPickle as pickle
import logging as logger

#refReadFile = "./forcredit/hw1_W_2/ref_hw1_W_2_chr_1.txt"
#readFile = "./forcredit/hw1_W_2/reads_hw1_W_2_chr_1.txt"

refReadFile = "../HW2/practice_E_1/ref_practice_E_1_chr_1.txt"
readFile = "../HW2/practice_E_1/reads_practice_E_1_chr_1.txt"

unmatchedReads = tools.readRead(readFile)

refSeq = tools.ReferenceSequence(refReadFile)
a = tools.generateKmerMap(refSeq.refRead,len(unmatchedReads[0].read))
i =0
len_unmatchedReads = len(unmatchedReads)
logger.warning("start sequencing")

for read in unmatchedReads:
    i+=1
    if i % (len_unmatchedReads/100) == 0:
        logger.warning("we have completed " + str(100*i/float(len_unmatchedReads)) + "% of this read")
    if refSeq.findMatch(read,a) == False:
        if refSeq.findMatch(read,a,reverse = True) == False:
            if read.pairedRead.matchRangeInRef is not None:
                #The pairedEnd had a match So lets check within a close distance
                start = read.pairedRead.matchRangeInRef[1]+50
                refSeq.findInDels(read,start,start + 400)
            else:
                #the pairedEnd did not find a match so lets check over the whole distance
                #refSeq.findInDels(read,0,len(refSeq.refRead))
                pass
          
refSeq.printInfo()
