'''
Created on Jan 8, 2016

@author: dulupinar
'''
import sequencingtools.tools as tools
import numpy as np
import os.path
import cPickle as pickle
import logging as logger
import operator
from operator import itemgetter
from _collections import defaultdict

def reject_outliers(data, m=2):
    return data[abs(data - np.mean(data)) < m * np.std(data)]

#refReadFile = "./forcredit/hw1_W_2/ref_hw1_W_2_chr_1.txt"
#readFile = "./forcredit/hw1_W_2/reads_hw1_W_2_chr_1.txt"
output = open("str_test_grad_write_cov.75.txt",'w') 
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
refSeq = tools.ReferenceSequence(refReadFile)
kmerMap = tools.generateKmerMap(refSeq.refRead,50)
    
#####
unmatchedReads = tools.readRead(readFile)
i =0
len_unmatchedReads = len(unmatchedReads)


if True:
    #this is the bad way to do everything
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
    
        
        #now we check for STR
        if (bool(found_read1) ^ bool(found_read2)):
            paired_reads = []
            if found_read1 == False:
                if read2_reverse == False:
                    paired_reads.append(read[0][::-1])
                else:
                    paired_reads.append(read[0])
                    
            if found_read2 == False:
                if read1_reverse == False:
                    paired_reads.append(read[1][::-1])
                else:
                    paired_reads.append(read[1])
                
            for read in paired_reads:
                for key in all_kmer:
                    if  all_kmer[key].matchKmer(read) == True:
                        if found_read1 == False :
                            starter = found_read2
                            if read2_reverse == True:
                                r = read[1][::-1]
                            else:
                                r = read[1]
                        else:
                            starter = found_read1
                            if read1_reverse == True:
                                r = read[0][::-1]
                            else:
                                r = read[0]
                        all_kmer[key].start_pos.append(starter)
                        all_kmer[key].matchKmer(r)


    for key in sorted(all_kmer):
        end_points = []
        size_count = []
        for read in all_kmer[key].reads:
            size_count.append(read.count(key))
            start_point = read.find(key*3)
            if start_point > 12:
                pos = kmerMap[read[start_point-12:start_point]]
                if len(pos) == 1:
                    end_points.append(pos)
        
            end_point = read.rfind(key*3)+len(key)*3
            if len(read) - end_point > 12:
                pos = kmerMap[read[end_point:end_point + 12]]
                if len(pos) == 1: 
                    end_points.append(pos)
                    
        if end_points != [] :
            #change this so we know the orientation so we get correct orienation
            #end_points = reject_outliers(np.array(end_points), m = 3)
            end_points = [end_point for end_point in end_points[0] if abs(end_point-all_kmer[key].start_pos[-1]) < 500]
            #filter_dict = defaultdict(list)
            #for point in end_points:
            #    filter_dict[str(point)[0:2]].append(point)        
            #max_l = -1
            #best_key = ""
            #summ = 0
            #for key_1 in filter_dict:
            #    if len(filter_dict[key_1]) > max_l:
            #        max_l = len(filter_dict[key_1])
            #        best_key = key_1
            #    summ+=len(filter_dict[key_1])
                    
            #end_points = np.array(filter_dict[best_key])
            end_points = np.asarray(end_points)
            std_dev = end_points.std()
            all_kmer[key].std_dev = std_dev
            all_kmer[key].end_points = end_points
            #all_kmer[key].avg_size = sum(size_count)/len(size_count)
            #all_kmer[key].filter_dict = filter_dict
    if True:               
        set_all = set()
        for key in all_kmer:
            if all_kmer[key].std_dev < 20 and all_kmer[key].std_dev != -1:
                set_all.add((all_kmer[key].end_points.min(),all_kmer[key].end_points.max()-all_kmer[key].end_points.min(),key))
                
        set_all_sorted = sorted(set_all,key=itemgetter(0))
        final_set = set()
        final_set.add(set_all_sorted[0])
        
        for i in range(len(set_all_sorted) - 1)[1:]:
            if set_all_sorted[i][0] > set_all_sorted[i -1][0] + 500:
                final_set.add(set_all_sorted[i])
        
        for s in final_set:
            refSeq.STR.add((s[0],s[2]*10))
            
refSeq.printInfo(output)

