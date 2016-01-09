'''
Created on Jan 8, 2016

@author: dulupinar
'''
import unittest
import sequencingtools.tools as tools
class TestSequencer1(unittest.TestCase):
    
    def test_del1(self):
        refReadFile = "./practice_W_1/ref_practice_W_1_chr_1.txt"
        readFile = "./practice_W_1/reads_practice_W_1_chr_1.txt"
        unmatchedReads = tools.readRead(readFile)
        refSeqObj = tools.ReferenceSequence(refReadFile)        
        '''
        for read in unmatchedReads:
            refSeqObj.findInDel(read,0,len(refSeqObj.refRead))
        refSeqObj.printInfo()
        '''
 
class TestSequencer2(unittest.TestCase):
           
    def test_del2(self):
        refRead = 'ACCAACCACGCCATGCTATAGGAGGGCCTTACGGGGACGAATGGGGAGACCAAACTTCTATGTTCGAGAGTACCGATATATCTCATTAAGAGCCGAGCAG'
        #refRead = 'CCATGCGCGGTGCGTATACCCCCATGCTAAAGCTAGACTACTAGTGATACAATCTTTTGCTGAGCTAGCGGCTTGTGCTACTAACTAAAAGTCCCAGGAT'
        readFile = "./practice_W_1/reads_practice_W_1_chr_1.txt"
        unmatchedReads = tools.readRead(readFile)
        refSeqObj = tools.ReferenceSequence(refRead,debug=True)
        
        for read in unmatchedReads:
            if refSeqObj.findInDel(read,0,len(refRead)):
                print unmatchedReads.index(read)
        refSeqObj.printInfo()

class TestSequencer3(unittest.TestCase):   
    def test_del3(self):
        
        refRead = "xxxxxTTAppppppppppppppppppppppppppppppppppppppppppppACCTTACCAAGTCCAACGATGATTGTCGCGCT"
        read    = "TTACCAAGTGATGATTGTCGCGCT"  # took out ACCTTACCAAGTCCAACGATGATTGTCGCGCT
        
        kmerMap = tools.generateKmerMap(refRead,len(read))
        refSeqObj = tools.ReferenceSequence(refRead,debug=True)
        readObj = tools.ReadSequence(read)
        refSeqObj.findInDel(readObj,3,len(refRead))
        refSeqObj.printInfo()

class TestSequencer4(unittest.TestCase):   
    def test_del(self):
        
        refReadFile = "./practice_W_1/ref_practice_W_1_chr_1.txt"
        readFile = "./practice_W_1/reads_practice_W_1_chr_1.txt"
        unmatchedReads = tools.readRead(readFile)
        refSeq = tools.ReferenceSequence(refReadFile)
        a = tools.generateKmerMap(refSeq.refRead,len(unmatchedReads[0].read))
        val = refSeq.findMatch(unmatchedReads[977],a,reverse=True)
        print val
        
if __name__ == "__main__":

    suite = unittest.TestLoader().loadTestsFromTestCase(TestSequencer4)
    unittest.TextTestRunner().run(suite)
