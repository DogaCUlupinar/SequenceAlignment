'''
Created on Jan 8, 2016

@author: dulupinar
'''
import sequencingtools.tools as tools
refReadFile = "./practice_W_1/ref_practice_W_1_chr_1.txt"
read = "CGTATTAGGAAAACGGTGTAAGGAGTAAAGCCGGTAGGGGAGTCTGACCG"
readFile = "./practice_W_1/reads_practice_W_1_chr_1.txt"
reads = tools.readRead(readFile)

refSeq = tools.ReferenceSequence(refReadFile)
a = tools.generateKmerMap(refSeq.refRead,len(read))

for read in reads:
    if refSeq.findMatch(read,a) == False:
        print "trying reverse"
        refSeq.findMatch(read[::-1],a)
    
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