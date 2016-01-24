#!/usr/local/bin/python
'''
Created on Jan 8, 2016

@author: dulupinar
'''
from multiprocessing import Process, Manager, cpu_count, Queue
import sequencingtools.tools as tools
import gc
import logging as logger
import glob
import sys
from __builtin__ import file

logger.basicConfig(level=logger.WARNING,format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
q = Queue()
manager = Manager()

#refReadFile = "../HW2/practice_E_1/ref_practice_E_1_chr_1.txt"
                    
refReadFile = "../HW2/practice_W_1/ref_practice_W_1_chr_1.txt"
readPrefix = "../HW2/practice_W_1/reads/splitreads*"

refReadFile = "../HW2/hw2grad/ref_hw2grad_M_1_chr_1.txt"
readPrefix = "../HW2/hw2grad/reads/splitreads*"


refSeq = tools.ReferenceSequenceSafe(refReadFile,manager)
logger.warning("reading in reference file: " + refReadFile)
a = tools.generateKmerMap(refSeq.refRead,50,manager)

#refReadFile = "../HW2/hw2grad/ref_hw2grad_M_1_chr_1.txt"
#readPrefix = "../HW2/hw2grad/reads/splitreads*"


readFiles = [file for file in glob.glob(readPrefix)]
for file in readFiles:
    q.put(file)
    


#create list of 2 workers
workers = [Process(target=tools.processRead, args=(x,q,a,refSeq)) for x in range(cpu_count())]
for worker in workers:
    worker.start()
for worker in workers:
    worker.join()

aligner = tools.Aligner(refSeq)
aligner.printInfo()
