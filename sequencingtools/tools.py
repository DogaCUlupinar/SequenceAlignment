'''
Created on Jan 7, 2016

@author: dulupinar
'''

DEFAULT_BLOCK_SIZE = 4
DEFAULT_ERRORS     = 1                 #errors per Read
DEFAULT_NUM_BLOCKS = DEFAULT_ERRORS +1 #kmer + 1
DEFAULT_COVERAGE = 5

class Variation:
    'SNP sites'
    def __init__(self,pos,variation):
        self.pos = pos
        self.var = variation
        
    def __repr__(self):
        return "Position: {pos},Variation: {var}".format(pos=self.pos,var=self.var)
    
    def __str__(self):
        return "Position: {pos},Variation: {var}".format(pos=self.pos,var=self.var)
        
class ReferenceSequence:
    'Reference sequence class'
    def __init__(self,filename):
        self.refRead = readRef(filename)
        self.variations = dict()
        self.coverage = []
        self.SNP = []
        
    def addVariations(self,readVariations):
        for variation in readVariations:
            if self.variations.has_key(variation.pos):
                self.variations[variation.pos].append(variation.var)
            else:
                self.variations[variation.pos] = [variation.var]
    
    def determineSNP(self,coverage = DEFAULT_COVERAGE):
        for key in self.variations:
            if len(self.variations[key]) > coverage:
                snp = max(self.variations[key],key=self.variations[key].count)
                self.SNP.append((self.refRead[key],snp,key))
                
    def findMatch(self,read,kmerMap,numBlocks=DEFAULT_NUM_BLOCKS):
        
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
                        return True
        
        return False
    
    def printInfo(self):
        self.determineSNP()
        print self.SNP
        
    def updateCoverage(self,start,end):
        pass

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
    return refRead

def readRead(refFilename):
    paired_end_reads = []
    with open(refFilename) as f:
        descriptorLine = f.readline()
        if ">" not in descriptorLine:
            raise Exception("Did not find > in first line are you sure this is read file")
        for line in f:
            #how do we want to handle paired ends
            paired_end_reads+=(line.strip().split(',')) # The two paired ends are separated by a comma
    return paired_end_reads


'''
def main():
    refRead = "ACCTTACCAAGTCCAACGATGATTGTCGCGCTGTCGGTAACGAGGAAAGAATCAGAATGCTAAGAGAATACCGAACCTAC"
    readLength = 12
    read = "TCCAAGGATGGT"
    numMerBlocks = 3
    
    a = generateKmerMap(refRead,len(read),numMerBlocks)
    refSeq = ReferenceSequence(refRead)
    refSeq.findMatch(read,a)
    refSeq.printInfo()

if __name__ == "__main__":
    main()
'''