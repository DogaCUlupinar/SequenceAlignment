'''
Created on Jan 7, 2016

@author: dulupinar
'''
from __future__ import print_function
from bisect import bisect_left
from collections import defaultdict, Counter, OrderedDict
import numpy as np
import cPickle as pickle
from sets import Set
from operator import itemgetter
from sys import stdout
import re
import logging as logger
import os
import sys
import operator

#from guppy import hpy
#hp = hpy()
alpha = ["A","C","T","G"]


logger.basicConfig(level=logger.WARNING,format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


DEFAULT_ERRORS     = 3                 #errors per Read
DEFAULT_NUM_BLOCKS = DEFAULT_ERRORS +1 #kmer + 1
DEFAULT_COVERAGE   = 8


class CNV(set):
    def __init__(self):
        self.header = ">CNV"
        set.__init__(self)
    
    def __str__(self):
        return_string = self.header + "\n"
        return return_string
    
class ALU(set):
    def __init__(self):
        self.header = ">ALU"
        set.__init__(self)
        
    def __str__(self):
        return_string = self.header + "\n"
        return return_string 
class INV(set):
    def __init__(self):
        self.header = ">INV"
        set.__init__(self)
        
    def __str__(self):
        return_string = self.header + "\n"
        return return_string
    
class STR(set):
    def __init__(self):
        self.header = ">STR"
        set.__init__(self)
        
    def __str__(self):
        return_string = self.header + "\n"
        for tup in sorted(self,key=itemgetter(0)):
            return_string+=tup[1]+","+str(tup[0])+"\n"
        return return_string
    
        
class SNP(set):
    def __init__(self):
        self.header = ">SNP"
        set.__init__(self)
    
    def __str__(self):
        return_string = self.header + "\n"
        for tup in sorted(self,key=itemgetter(2)):
            return_string+=tup[0] + "," + tup[1] + "," + str(tup[2]) + "\n"
        return return_string

class INS(set):
    def __init__(self):
        self.header = ">INS"
        set.__init__(self)
    
    def __str__(self):
        return_string = self.header + "\n"
        for tup in sorted(self,key=itemgetter(1)):
            return_string+=tup[0] + "," + str(tup[1]) + "\n"
        return return_string
    
class DEL(set):
    def __init__(self):
        self.header = ">DEL"
        set.__init__(self)
    
    def __str__(self):
        return_string = self.header + "\n"
        for tup in sorted(self,key=itemgetter(1)):
            return_string+=tup[0] + "," + str(tup[1]) + "\n"
        return return_string
    
class Variation:
    'SNP sites'
    def __init__(self,pos,variation):
        self.pos = pos
        self.var = variation
        
    def __repr__(self):
        return "Position: {pos},Variation: {var}".format(pos=self.pos,var=self.var)
    
    def __str__(self):
        return "Position: {pos},Variation: {var}".format(pos=self.pos,var=self.var)
    
class ReadSequence:
    'Read sites has information like where in the reference and its paired end'
    
    def __init__(self,read): 
        self.read = read
        self.pairedRead = None
        self.matchRangeInRef = None
        
    def setMatchRangeInRef(self,start,end):
        self.matchRangeInRef = (start,end)
        
    def __repr__(self):
        return self.read
    
    def __len__(self):
        return len(self.read)
        
class ReferenceSequence:
    'Reference sequence class'
    def __init__(self,filename,debug=False):
        if debug:
            #the filename just becomes the refRead makes debugging easier
            self.refRead = filename 
            self.name = "DEBUG"
        else:
            fileContent = readRef(filename)
            self.refRead = fileContent[0]
            self.name = fileContent[1]
        
        #attributes
        self.STR = STR()
        self.CNV = CNV()
        self.ALU = ALU()
        self.INV = INV()
        self.insertions = INS()
        self.deletions = DEL()
        self.SNP = SNP()
        self.attributes = [self.name,self.STR,self.CNV,self.ALU,self.INV,self.insertions,self.deletions,self.SNP]
        
        #helper
        self.variations = dict()
    def generateDonor(self):
        self.determineSNP()
        self.determineInDels()
        self.donor_seq = list(self.refRead)
        self.index_match  = np.linspace(0,len(self.refRead) -1,num=len(self.refRead),dtype=int) #key reference position value is position in ref
        for snp in self.SNP:
            self.donor_seq[snp[2]] = snp[1]
        
        index_offset = 0
        match_offset = 0
        last_update = 0
        #insertions are find but when i use them my score goes to poop
        for ins in sorted(self.insertions,key=itemgetter(1)):
            #update index_match
            for i in range(last_update,ins[1]):
                self.index_match[i]+=match_offset
                
            self.donor_seq.insert(ins[1] + index_offset,ins[0])
            index_offset += 1
            match_offset+=len(ins[0])
            last_update = ins[1]
        
        for i in range(last_update,len(self.refRead)):
                self.index_match[i]+=index_offset
        
        """
        for dell in self.deletions:
            del self.donor_seq[dell[1] + index_offset:dell[1]+len(dell[0]) + index_offset]
            index_offset -=len(dell[0])
        Over it
        """
        
        self.donor_seq = "".join(self.donor_seq)
        return self.donor_seq
    
    def translateIndex(self,donor_index):
        while True:
            try:
                match = np.where(self.index_match==donor_index)[0][0]
                break
            except IndexError:
                donor_index+=1
        
        return match      
    def findSTRRegex(self):
        all_kmer = dict()
        allKmer(3, "", alpha, all_kmer)
        allKmer(4, "", alpha, all_kmer)
        allKmer(5, "", alpha, all_kmer)
        
        len_dict = len(all_kmer)
        ii = 0
        dict_str = dict() #key is start pos, #val is string
        for key in all_kmer:
            if ii % (len_dict/20) == 0:logger.warn("Checking the {0}th kmer".format(str(ii)))
            ii+=1
            for m in re.finditer("({kmer}){{{rep},}}".format(kmer=key,rep=4-(len(key)-4)),self.donor_seq):
                match_inref = self.translateIndex(m.start())
                dict_str[match_inref] = m.group()
               
        sorted_d = sorted(dict_str.items(), key=operator.itemgetter(0))
        self.STR.add(sorted_d[0])
        
        for index in range(len(sorted_d) -1)[1:]:
            prev = sorted_d[index -1][0]
            current = sorted_d[index][0]
            if prev < current - 5:
                self.STR.add(sorted_d[index])
                
        return self.STR
        
    def findSTR(self):
        #find all str
        str_count = defaultdict(list) #key is kmer value is postions
        str_dict  = dict() #key is position value is str
        repeat_min = 4
        len_read = len(self.refRead)
        for i in range(len(self.refRead) - 6):
            if i % (len_read/25) == 0:logger.warn("processed {0} of {1} of read".format(i,len_read))
            kmer_3 = self.refRead[i:i+3]
            kmer_4 = self.refRead[i:i+4]
            kmer_5 = self.refRead[i:i+5]
            str_count[kmer_3].append(i)
            str_count[kmer_4].append(i)
            str_count[kmer_5].append(i)
        
        ii = 0
        len_dict = len(str_count)
        for kmer in str_count:
            if ii % (len_dict/20) == 0:logger.warn("Checking the {0}th kmer".format(str(ii)))
            ii+=1         
            len_kmer = len(kmer)
            positions = str_count[kmer]
            
            #for each str
            pos_str = -1 #position of str
            str_seq = "" #str itself
            repeat_length = 1 #number of times str is repeated
            index = 0
            while index in range(len(positions) -1):
                pos = positions[index]
                find = getIndex(positions,pos+len_kmer,index)
                if find != False:
                    #finding str
                    if pos_str == -1:
                        pos_str = pos                  
                    repeat_length+=1
                    str_seq+=kmer
                    index = find
                else:
                    if repeat_length > repeat_min:
                        str_dict[pos_str] = str_seq
                    pos_str = -1
                    str_seq = ""
                    repeat_length = 1
                    index+=1
                
            if repeat_length > repeat_min:
                    str_dict[pos_str] = str_seq
            
        sorted_d = sorted(str_dict.items(), key=operator.itemgetter(0))
        self.STR.add(sorted_d[0])
        
        for index in range(len(sorted_d) -1)[1:]:
            prev = sorted_d[index -1][0]
            current = sorted_d[index][0]
            if prev < current - 5:
                self.STR.add(sorted_d[index])
         
        return str_dict,str_count     
            
        
    def addVariations(self,readVariations):
        for variation in readVariations:
            if self.variations.has_key(variation.pos):
                list = self.variations[variation.pos]
                list.append(variation.var)
                self.variations[variation.pos] = list
            else:
                self.variations[variation.pos] = [variation.var]
    
    def determineSNP(self,coverage = DEFAULT_COVERAGE):
        for key in self.variations:
            snp = max(self.variations[key],key=self.variations[key].count)
            if self.variations[key].count(snp) >= coverage:
                self.SNP.add((self.refRead[key],snp,key))
                
    def findMatch(self,read,kmerMap,numBlocks=DEFAULT_NUM_BLOCKS,reverse=False):
        if reverse:
            read = read[::-1]
        blockSize = len(read)/numBlocks
        for explore_i in range(numBlocks):
            #check if any block matches our preproccessed map
            start = explore_i*blockSize
            end = min(len(read),start+blockSize) #make better need to get read of len and min operators
            block = read[start:end]
            if kmerMap.has_key(block):
                for posInRef in kmerMap[block]:
                    #verify around each site where the block matches
                    startRef = -explore_i*blockSize + posInRef
                    endRef   = min(len(self.refRead),startRef+len(read)) #make better need to get read of len and min operators
                    variations = diffPlaces(self.refRead,read,startRef,endRef)
                    if len(variations) <= numBlocks-1:
                        self.addVariations(variations)
                        #self.updateCoverage(startRef,endRef)
                        return startRef
        
        return False
    
    def findInDels(self,readObj,start_ref,end_ref):
        if self.findInsertions(readObj,start_ref,end_ref) == False:      
            self.findDeletions(readObj,start_ref,end_ref)
        
    def findDeletions(self,read,start_ref,end_ref):
        refRead = self.refRead[start_ref:end_ref]
        'currently for deletions'
        last_match_index = 0

        while last_match_index <= len(refRead)-len(read):
            logger.info("last match: {0}, length of refRead {1}, length of read {2}".format(str(last_match_index),str(len(refRead)),str(len(read))))
            deletion = ""
            i = 0
            j = last_match_index
            last_match_index+=1
            matches = 0
            gap_size = 0
            gap_count = 0
            max_gap_count = 1
            max_gap_size = 20
            match_previous = False
            gap_pos = None
            
            while gap_size <= max_gap_size and gap_count <= max_gap_count and i < len(read) and j < len(refRead):
                logger.info("gap size {0}, gap count {1}, i {2}, j{3}".format(str(gap_size), str(gap_count), str(i), str(j)))
                if refRead[j] == read[i]:
                    j+=1
                    i+=1
                    matches+=1
                    match_previous = True
                    last_match_index = j
                    
                else:
                    if match_previous:
                        if refRead[j] == read[i-1]:
                            match_previous = True
                        else:
                            gap_count+=1
                            gap_pos = j
                            match_previous = False
                        
                    if gap_count > 0:
                        #we are in a gap
                        #only add as deletion after we have started matching the read to reference
                        gap_size+=1
                        deletion+=refRead[j]
                    j+=1 #move j pointer but i stays same

            
            if matches == len(read) and gap_count > 0 and gap_count <=max_gap_count:
                check = gap_pos + gap_size -1
                if len(read) == last_match_index - check:
                    #this is ugly but need to skip these
                    continue
                start = start_ref+gap_pos
                self.updateDeletions(start,deletion)
                #self.updateCoverage(start_ref, end_ref)
                return True

        return False
    
    def findInsertions(self,read,start_ref,end_ref):
        refRead = self.refRead[start_ref:end_ref]
        'currently for Insertions'
        last_match_index = 0

        while last_match_index <= len(refRead)-len(read):
            insertion = ""
            i = 0
            j = last_match_index
            matches = 0
            gap_size = 0
            gap_count = 0
            max_gap_count = 1
            max_gap_size = 20
            match_previous = False
            gap_pos = None
            last_match_index+=1
            
            while gap_size <= max_gap_size and gap_count <= max_gap_count and i < len(read) and j < len(refRead):
                if refRead[j] == read[i]:
                    j+=1
                    i+=1
                    matches+=1
                    match_previous = True
                    last_match_index = j
                    
                else:
                    if match_previous:
                        if refRead[j] == read[i-1]:
                            match_previous = True
                        else:
                            gap_count+=1
                            gap_pos = j
                            match_previous = False
                        
                    if gap_count > 0:
                        #we are in a gap
                        #only add as deletion after we have started matching the read to reference
                        gap_size+=1
                        insertion+=read[i]
                    i+=1 #move j pointer but i stays same

            
            if matches == len(read) - gap_size and gap_count > 0 and gap_count <=max_gap_count:
                check = gap_pos + gap_size -1
                if len(read) == last_match_index - check:
                    #this is ugly but need to skip these
                    continue
                start = start_ref+gap_pos
                self.updateInsertions(start,insertion)
                return True

        return False
    
    def printInfo(self,filestream=stdout):
        self.determineSNP()
        self.determineInDels()
        for attribute in self.attributes:
            filestream.write(str(attribute))

        
    def updateCoverage(self,start,end):
        for i in range(start,end):
            self.coverage[i]+=1

    def updateDeletions(self,start,deletion):
        self.deletions.add((deletion,start))
            
    def updateInsertions(self,start,insertion):
        self.insertions.add((insertion,start))
    
    def determineInDels(self):
        self.INDEL = self.deletions | self.insertions
        remove = Set()
        for dely in self.deletions:
            for insy in self.insertions:
                dist = abs(dely[1] - insy[1])
                max_len = max(len(dely[0]),len(insy[0]))
                if dist < max_len:
                    if len(dely[0])>len(insy[0]):

                        remove.add(dely)
                    else:
                        remove.add(insy)
 
        for tup in remove:
            if tup in self.deletions:
                print("WARNING we are removing dely " + str(tup))
                self.deletions.remove(tup)
            if tup in self.insertions:
                print("WARNING we are removing insy " + str(tup))
                self.insertions.remove(tup)

               
    def determineDeletions(self):
        deletion = ""
        
        for pos,value in enumerate(self.deletions):
            
            if value == 0:
                if deletion != "":
                    self.INDEL.add((deletion,pos-len(deletion)))
                    deletion = ""
            else:         
                deletion += self.refRead[pos]

def getIndex(a, x,start):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(a, x,lo=start)
    if i != len(a) and a[i] == x:
        return i
    else:
        return False
                   
def diffPlaces(refRead,read,start,end):
    #check the places where str1 and str2 are different
    #currently string should be of the same length
    numDifferences = []
    blockRef = refRead[start:end]
    for i in range(len(blockRef)):
        if blockRef[i] != read[i]:
            #pos,original,diff
            numDifferences.append(Variation(start+i,read[i]))
            
    return numDifferences
   
def generateKmerMap(refRead, readLength, numMerBlocks=DEFAULT_NUM_BLOCKS):
    '''
    @param numMerBlocks is the number chunks we want to divide the reads into; kmer + 1
    '''
    refReadLengh = len(refRead)
        
    pickle_file = "kmerMap_{0}_{1}.p".format(numMerBlocks,refReadLengh)
    if os.path.isfile(pickle_file):
        logger.warning("reading in KmerMap from " + pickle_file)
        return pickle.load(open(pickle_file,"rb"))
    else:
        logger.warning("generating Map")
        chunkSize = min(readLength,readLength/numMerBlocks) #the size of each kmer block
        #hp.setrelheap()
        kmerMap = defaultdict(list)
        rangeR = refReadLengh-chunkSize

        for i in range(rangeR):
            if i % (rangeR/15) == 0:
                logger.warning("generated {0}% of the kmerMap with size {1}".format(str(100*i/float(rangeR)),str(sys.getsizeof(kmerMap)))) 
                #h = hp.heap()
                #print(h)
            start = i
            end = start + chunkSize #just ignore if it is not correct size
            block = refRead[start:end]
            kmerMap[block].append(start)
                
        logger.warning("creating pickle file {0}".format(pickle_file))
        #pickle.dump(kmerMap,open(pickle_file,"wb"))
        return kmerMap
    
class STRCandidate():
    def __init__(self,pattern):
        self.pattern = pattern #this is read
        self.reads = []
        self.std_dev = -1
    
    def matchKmerBAD(self,read):
        
        match = self.pattern.match(read)
        if match:
            self.addRead(read)
            return True
        else:
            return False
        
    def matchKmer(self,read):
             
        if self.pattern*5 in read:
            self.addRead(read)
            return True
        else:
            return False
        
    def addRead(self,read):
        self.reads.append(read)
    
def allKmerBetter(end,kmer,alpha,counter):
    if len(kmer) == end:
        #pattern = re.compile("\w{{0,24}}((?:{kmer_pat}){{3,20}})\w{{0,24}}".format(kmer_pat=kmer))
        
        counter[kmer]= STRCandidate(kmer)
        return counter

    for i in alpha:
        kmer_new= kmer + i
        allKmerBetter(end,kmer_new,alpha,counter)
    
    return counter

def allKmer(end,kmer,alpha,counter):
    if len(kmer) == end:
        #pattern = re.compile("\w{{0,24}}((?:{kmer_pat}){{3,20}})\w{{0,24}}".format(kmer_pat=kmer))
        
        counter[kmer]=[]
        return counter

    for i in alpha:
        kmer_new= kmer + i
        allKmer(end,kmer_new,alpha,counter)
    
    return counter
                
def stringDnaToNum(stringDna):
    number = ""
    for letter in stringDna:
        if letter == "A":
            number+="0"
        elif letter == "C":
            number+="1"
        elif letter == "T":
            number+="2"
        else:
            number+="3"
    return int(number)
    
def readRef(refFilename):
    refRead = ""
    with open(refFilename) as f:
        descriptorLine = f.readline()
        if ">" not in descriptorLine:
            raise Exception("Did not find > in first line are you sure this is reference file")
        for line in f:
            refRead += line.strip()
    return refRead,descriptorLine

def readRead(refFilename):
    paired_end_reads = []
    with open(refFilename) as f:
        descriptorLine = f.readline()
        if ">" not in descriptorLine:
            paired_end_read = descriptorLine.strip().split(',') # The two paired ends are separated by a comma

            if paired_end_read[0] != ""  and paired_end_read[1] != "":
                paired_end_reads.append((paired_end_read[0],paired_end_read[1]))
                
        for line in f:
            #how do we want to handle paired ends
            paired_end_read = line.strip().split(',') # The two paired ends are separated by a comma
            if len(paired_end_read[0]) < 40  or len(paired_end_read[1]) < 40:
                continue
            paired_end_reads.append((paired_end_read[0],paired_end_read[1]))
    
    return paired_end_reads

def processRead(id,queue,kmerMap,refSeq):
    read_count = 0
    while queue.empty() == False:
        readFile = queue.get()
        logger.warning("Process {0} reading in read file: {1}".format(str(id),readFile))
        unmatchedReads = readRead(readFile)    
        len_unmatchedReads = len(unmatchedReads)
       
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
'''        
def rr(all_kmer,kmer,corr):
    summ = 0
    for r in all_kmer[kmer].reads:
        summ+=r.count(kmer)
        
    print "summ of kmer: ", summ
    print "times in corr" ,corr.count(kmer)
    print "avg times", summ/corr.count(kmer)
'''

if __name__ == "__main__":
    dicter = defaultdict(list)
    a = allKmerBetter(3, "", alpha, dicter)
    print(a)