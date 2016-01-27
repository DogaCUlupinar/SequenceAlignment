#!/usr/local/bin/python
'''
Created on Jan 8, 2016

@author: dulupinar
'''
from multiprocessing import Process, Manager, cpu_count, Queue
import sequencingtools.tools as tools
import logging as logger
import glob


logger.basicConfig(level=logger.WARNING,format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
q = Queue()
manager = Manager()

#refReadFile = "../HW2/practice_E_1/ref_practice_E_1_chr_1.txt"
#readPrefix = "../HW2/practice_E_1/reads/splitreads*"
                    
#refReadFile = "../HW2/practice_W_1/ref_practice_W_1_chr_1.txt"
#readPrefix = "../HW2/practice_W_1/reads/splitreads*"

refReadFile = "../HW2lnx/hw2grad/ref_hw2grad_M_1_chr_1.txt"
readPrefix = "../HW2lnx/hw2grad/reads/splitreads*"
output = open("hw2practice.txt",'w')

refSeq = tools.ReferenceSequenceSafe(refReadFile,manager)
logger.warning("reading in reference file: " + refReadFile)
a = tools.generateKmerMap(refSeq.refRead,50,manager)

#refReadFile = "../HW2/hw2grad/ref_hw2grad_M_1_chr_1.txt"
#readPrefix = "../HW2/hw2grad/reads/splitreads*"


readFiles = [readfile for readfile in glob.glob(readPrefix)]
for readfile in readFiles:
    q.put(readfile)
    


#create list of 2 workers
workers = [Process(target=tools.processRead, args=(x,q,a,refSeq)) for x in range(cpu_count())]
for worker in workers:
    worker.start()
for worker in workers:
    worker.join()

aligner = tools.Aligner(refSeq)
aligner.printInfo(filestream=output)
