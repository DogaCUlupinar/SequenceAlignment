'''
Created on Jan 7, 2016

@author: dulupinar
'''
from __future__ import print_function
import cPickle as pickle
from sets import Set
from operator import itemgetter
from sys import stdout
import logging as logger
import os
import sys
#from guppy import hpy
#hp = hpy()

logger.basicConfig(level=logger.WARNING,format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


DEFAULT_ERRORS     = 3                 #errors per Read
DEFAULT_NUM_BLOCKS = DEFAULT_ERRORS +1 #kmer + 1
DEFAULT_COVERAGE   = 8

class STR(set):
    def __init__(self):
        self.header = ">STR"
        set.__init__(self)
        
    def __str__(self):
        return_string = self.header + "\n"
        return return_string
    
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
        
class Aligner:
    def __init__(self,refSeq):
        self.STR = STR()
        self.CNV = CNV()
        self.ALU = ALU()
        self.INV = INV()
        self.INS = INS()
        self.DEL = DEL()
        self.SNP = SNP()
        self.attributes = [self.STR,self.CNV,self.ALU,self.INV,self.INS,self.DEL,self.SNP]
        
        for item in refSeq.insertions:
            self.INS.add(item)
        for item in refSeq.deletions:
            self.DEL.add(item)
        
        #helper         
        self.determineInDels()
        self.determineSNP(refSeq.refRead,refSeq.variations)
        
    def determineInDels(self):

        remove = Set()
        for dely in self.DEL:
            for insy in self.INS:
                dist = abs(dely[1] - insy[1])
                max_len = max(len(dely[0]),len(insy[0]))
                if dist < max_len:
                    if len(dely[0])>len(insy[0]):

                        remove.add(dely)
                    else:
                        remove.add(insy)
 
        for tup in remove:
            if tup in self.DEL:
                print("WARNING we are removing dely " + str(tup))
                self.DEL.remove(tup)
            if tup in self.INS:
                print("WARNING we are removing insy " + str(tup))
                self.INS.remove(tup)
    
    def determineSNP(self,refRead,variations,coverage = DEFAULT_COVERAGE):
        var = dict(variations)
        for key in var:
            snp = max(variations[key],key=variations[key].count)
            if variations[key].count(snp) >= coverage:
                self.SNP.add((refRead[key],snp,key))
    
    def printInfo(self,filestream=stdout):
        for attribute in self.attributes:
            print(attribute,file=filestream,end="")
                
class ReferenceSequenceSafe:
    'Reference sequence class'
    def __init__(self,filename,manager,debug=False):
        self.manager = manager
                
        if debug:
            #the filename just becomes the refRead makes debugging easier
            self.refRead = filename 
        else:
            fileContent = readRef(filename)
            self.refRead = fileContent[0]
            self.name = fileContent[1]
        
        #attributes
        
        #helper
        self.insertions = manager.list()
        self.deletions = manager.list()

        self.variations = manager.dict()
        #self.coverage = np.zeros(len(self.refRead),dtype=np.int)

        
    def addVariations(self,readVariations):
        for variation in readVariations:
            if self.variations.has_key(variation.pos):
                list = self.variations[variation.pos]
                list.append(variation.var)
                self.variations[variation.pos] = list
            else:
                self.variations[variation.pos] = [variation.var]
                
    def findMatch(self,read,kmerMap,numBlocks=DEFAULT_NUM_BLOCKS,reverse=False):
    
        if reverse:
            read = read[::-1]
        blockSize = len(read)/numBlocks
        for explore_i in range(numBlocks):
            #check if any block matches our preproccessed map
            start = explore_i*blockSize
            end = min(len(read),start+blockSize) #make better need to get read of len and min operators
            block = read[start:end]
            block_num = stringDnaToNum(block)
            if kmerMap.has_key(block_num):
                for posInRef in kmerMap[block_num]:
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
    
    def updateDeletions(self,start,deletion):
        self.deletions.append((deletion,start))
            
    def updateInsertions(self,start,insertion):
        self.insertions.append((insertion,start))

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
   
def generateKmerMap(refRead, readLength,manager, numMerBlocks=DEFAULT_NUM_BLOCKS):
    
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
        kmerMap = dict()
        rangeR = refReadLengh-chunkSize

        for i in range(rangeR):
            if i % (rangeR/15) == 0:
                logger.warning("generated {0}% of the kmerMap with size {1}".format(str(100*i/float(rangeR)),str(sys.getsizeof(kmerMap)))) 
                #h = hp.heap()
                #print(h)
            start = i
            end = start + chunkSize #just ignore if it is not correct size
            
            block = refRead[start:end]
            block_num = stringDnaToNum(block)
            if kmerMap.has_key(block_num):
                kmerMap[block_num].append(start)
            else:
                kmerMap[block_num] = [start]
                
        #logger.warning("creating pickle file {0}".format(pickle_file))
        pickled_dict = pickle.dumps(kmerMap)
        print("THE FINAL SIZE" + str(sys.getsizeof(pickled_dict)))
        logger.warning("making the dictionary thread safe")
        kmerMap_safe = manager.dict(kmerMap)
        return kmerMap_safe
                
    
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
            paired_end_reads.append((paired_end_read[0],paired_end_read[1]))
        for line in f:
            #how do we want to handle paired ends
            paired_end_read = line.strip().split(',') # The two paired ends are separated by a comma
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
