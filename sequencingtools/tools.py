'''
Created on Jan 7, 2016

@author: dulupinar
'''
from __future__ import print_function
import numpy as np
from sets import Set
from operator import itemgetter
from sys import stdout


DEFAULT_ERRORS     = 2                 #errors per Read
DEFAULT_NUM_BLOCKS = DEFAULT_ERRORS +1 #kmer + 1
DEFAULT_COVERAGE   = 5

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
        
class ReferenceSequence:
    'Reference sequence class'
    def __init__(self,filename,debug=False):
        if debug:
            #the filename just becomes the refRead makes debugging easier
            self.refRead = filename 
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
        self.coverage = np.zeros(len(self.refRead),dtype=np.int)

        
    def addVariations(self,readVariations):
        for variation in readVariations:
            if self.variations.has_key(variation.pos):
                self.variations[variation.pos].append(variation.var)
            else:
                self.variations[variation.pos] = [variation.var]
    
    def determineSNP(self,coverage = DEFAULT_COVERAGE):
        for key in self.variations:
            snp = max(self.variations[key],key=self.variations[key].count)
            if self.variations[key].count(snp) >= coverage:
                self.SNP.add((self.refRead[key],snp,key))
                
    def findMatch(self,readObj,kmerMap,numBlocks=DEFAULT_NUM_BLOCKS,reverse=False):
        read = readObj.read
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
                    endRef   = min(self.refRead,startRef+len(read)) #make better need to get read of len and min operators
                    variations = diffPlaces(self.refRead,read,startRef,endRef)
                    if len(variations) <= numBlocks-1:
                        self.addVariations(variations)
                        self.updateCoverage(startRef,endRef)
                        readObj.setMatchRangeInRef(startRef,endRef)
                        return True
        
        return False
    
    def findInDels(self,readObj,start_ref,end_ref):
        self.findInsertions(readObj,start_ref,end_ref)
        self.findDeletions(readObj,start_ref,end_ref)
        
    def findDeletions(self,readObj,start_ref,end_ref):
        read = readObj.read
        refRead = self.refRead[start_ref:end_ref]
        'currently for deletions'
        last_match_index = 0

        while last_match_index <= len(refRead)-len(read):
            deletion = ""
            i = 0
            j = last_match_index
            matches = 0
            gap_size = 0
            gap_count = 0
            max_gap_count = 1
            max_gap_size = 20
            match_previous = False
            gap_pos = None
            
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
    
    def findInsertions(self,readObj,start_ref,end_ref):
        read = readObj.read
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
            print(attribute,file=filestream,end="")

        
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
    chunkSize = min(readLength,readLength/numMerBlocks) #the size of each kmer block
    kmerMap = dict()
    for i in range(refReadLengh-chunkSize):
        start = i
        end = min(len(refRead),start + chunkSize) #just ignore if it is not correct size
        block = refRead[start:end]
        if kmerMap.has_key(block):
            kmerMap[block].append(start)
        else:
            kmerMap[block] = [start]
    
    return kmerMap
                
    
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
            raise Exception("Did not find > in first line are you sure this is read file")
        for line in f:
            #how do we want to handle paired ends
            paired_end_read = line.strip().split(',') # The two paired ends are separated by a comma
            seq1 = ReadSequence(paired_end_read[0])
            seq2 = ReadSequence(paired_end_read[1])
            seq1.pairedRead = seq2
            seq2.pairedRead = seq1
            paired_end_reads.append(seq1)
            paired_end_reads.append(seq2)
            
    return paired_end_reads


