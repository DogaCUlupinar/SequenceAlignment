'''
Created on Jan 8, 2016

@author: dulupinar
'''
import sequencingtools.tools as tools

refReadFile = "./practice_W_1/ref_practice_W_1_chr_1.txt"
read = "CGTATTAGGAAAACGGTGTAAGGAGTAAAGCCGGTAGGGGAGTCTGACCG"
readFile = "./practice_W_1/reads_practice_W_1_chr_1.txt"
unmatchedReads = tools.readRead(readFile)

refSeq = tools.ReferenceSequence(refReadFile)
a = tools.generateKmerMap(refSeq.refRead,len(read))

for read in unmatchedReads:
    if refSeq.findMatch(read,a) == False:
        #print "trying reverse"
        if refSeq.findMatch(read,a,reverse = True):
            #print "found with reverse"
            pass
        else:
            #print "did not find chita choot"
            if read.pairedRead.matchRangeInRef is not None:
                start = read.pairedRead.matchRangeInRef[1]+50
                if refSeq.findInDel(read,start,start + 400):
                    print "found deletion"
                    
                else:
                    end = read.pairedRead.matchRangeInRef[1]
                    refSeq.findInDel(read,end - 400,end)
refSeq.printInfo()

'''
def main():
    refReadFile = "./practice_W_1/ref_practice_W_1_chr_1.txt"
    read = "CGTATTAGGAAAACGGTGTAAGGAGTAAAGCCGGTAGGGGAGTCTGACCG"
    
    refSeq = tools.ReferenceSequence(refReadFile)
    a = tools.generateKmerMap(refSeq.refRead,len(read))
    
    refSeq.findMatch(read,a)
    refSeq.printInfo()

if __name__ == "__main__":
    main()
'''