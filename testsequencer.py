'''
Created on Jan 8, 2016

@author: dulupinar
'''
import matplotlib.pyplot as plt
import unittest
import difflib
import sequencingtools.tools as tools
import sys
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
        refRead = 'TTACTTTCTGCGCGTCGGCTGAAGAGAGTGAGGACATTGATCTCCGCTGGATGGTCCCGGTCATGTAAACGCCGTTCTGGAGCCAATCCGTGCTCAGGGT'
        readFile = "./practice_W_1/reads_practice_W_1_chr_1.txt"
        unmatchedReads = tools.readRead(readFile)
        refSeqObj = tools.ReferenceSequence(refRead,debug=True)
        
        
        for read in unmatchedReads:
            a = refSeqObj.findInDel(read,0,len(refRead))
            if a is True:
                print unmatchedReads.index(read )
                
                
        refSeqObj.printInfo()

class TestSequencer9(unittest.TestCase):
           
    def test_del3(self):
        
        read       = "ACCTTACCAAGTAAAACCAACGATGATTGTCGCGCT"
        refRead    = "ACCTTACCAAGTCCAACGATGATTGTCGCGCTGGGGGGGGGGGGGGGG"  # took out ACCTTACCAAGTCCAACGATGATTGTCGCGCT
        
        kmerMap = tools.generateKmerMap(refRead,len(read))
        refSeqObj = tools.ReferenceSequence(refRead,debug=True)
        readObj = tools.ReadSequence(read)
        refSeqObj.findInDel2(readObj,0,len(refRead))
        refSeqObj.printInfo()
        
class TestSequencer5(unittest.TestCase):
           
    def test_del2(self):
        refRead = 'ACCAACCACGCCATGCTATAGGAGGGCCTTACGGGGACGAATGGGGAGACCAAACTTCTATGTTCGAGAGTACCGATATATCTCATTAAGAGCCGAGCAG'
        #refRead = 'CCATGCGCGGTGCGTATACCCCCATGCTAAAGCTAGACTACTAGTGATACAATCTTTTGCTGAGCTAGCGGCTTGTGCTACTAACTAAAAGTCCCAGGAT'
        refRead = 'GTCAGCTCTATCAGTATTGGTTTGTAGTGAAGAGAGATGTTTCTGCAATACCCCTGTTACTTACCCCTGTGTAGCAGGCGGCGTGAGCGAGAGCCATACAGGGTCCTTGATGCCGGAAACGGGTTACGTGCGAAGCTCCCGCAGCCGGTCCACTATCGTAAGGCCTTCCGTCGCAAGTAGTCATATTGCACGTTTTCTCG'
        readFile = "./practice_W_1/reads_practice_W_1_chr_1.txt"
        unmatchedReads = tools.readRead(readFile)
        refSeqObj = tools.ReferenceSequence(refRead,debug=True)
        
        refSeqObj.findInDel(unmatchedReads[1995],0,len(refRead))
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

class TestSequencer6(unittest.TestCase):   
    def test_del3(self):
        
        read =                                                         'TGGGTCCCGGTCATGTAAACGCCGTTCTGGAGCCAATCCGTGCTCAGGGT'
        refsnippet = 'TTACTTTCTGCGCGTCGGCTGAAGAGAGTGAGGACATTGATCTCCGCTGGATGGTCCCGGTCATGTAAACGCCGTTCTGGAGCCAATCCGTGCTCAGGGT'
        
        for i in range(len(read)):
            if read[i] != refsnippet[i]:
                print i
                print read[i]
                print refsnippet[i]
        
class TestSequencer4(unittest.TestCase):   
    def test_del(self):
        
        refReadFile = "./practice_W_1/ref_practice_W_1_chr_1.txt"
        readFile = "./practice_W_1/reads_practice_W_1_chr_1.txt"
        unmatchedReads = tools.readRead(readFile)
        refSeq = tools.ReferenceSequence(refReadFile)
        a = tools.generateKmerMap(refSeq.refRead,len(unmatchedReads[0].read))
        val = refSeq.findInDel(unmatchedReads[977],0,len(refSeq.refRead))
        print val
        refSeq.printInfo()
  
class TestSequencer8(unittest.TestCase):   
    def test_del(self):
        
        refReadFile = "./practice_W_1/ref_practice_W_1_chr_1.txt"
        readFile = "./practice_W_1/reads_practice_W_1_chr_1.txt"
        unmatchedReads = tools.readRead(readFile)
        refSeq = tools.ReferenceSequence(refReadFile)
        a = tools.generateKmerMap(refSeq.refRead,len(unmatchedReads[0].read))
        val = refSeq.findInDel(unmatchedReads[5220],0,len(refSeq.refRead))
        print val
        refSeq.printInfo()
        
class TestSequencerRegress(unittest.TestCase):   
    def test_del(self):
            refReadFile = "../HW1/practice_W_1/ref_practice_W_1_chr_1.txt"
            readFile = "../HW1/practice_W_1/reads_practice_W_1_chr_1.txt"
            gen_file = "../HW1/forcredit/testcase_gen.txt"
            solution_file = "../HW1/forcredit/testcase.txt"
            output = open(gen_file,'w')  
            unmatchedReads = tools.readRead(readFile)
            
            refSeq = tools.ReferenceSequence(refReadFile)
            a = tools.generateKmerMap(refSeq.refRead,len(unmatchedReads[0].read))
            
            for read in unmatchedReads:
                if refSeq.findMatch(read,a) == False:
                    if refSeq.findMatch(read,a,reverse = True) == False:
                        if read.pairedRead.matchRangeInRef is not None:
                            #The pairedEnd had a match So lets check within a close distance
                            start = read.pairedRead.matchRangeInRef[1]+50
                            refSeq.findInDels(read,start,start + 400)
                        else:
                            #the pairedEnd did not find a match so lets check over the whole distance
                            #refSeq.findInDels(read,0,len(refSeq.refRead))
                            pass
                      
            refSeq.printInfo(filestream=output)
            output.close()
            
            with open(solution_file,'r') as sol:
                with open(gen_file,'r') as gen:
                    diff = difflib.unified_diff(sol.readlines(), gen.readlines(), fromfile=solution_file, tofile=gen_file)
            
            Fail = False
            for line in diff:
                Fail = True
                sys.stdout.write(line)
                
            self.assertEqual(Fail,False)
                
if __name__ == "__main__":

    suite = unittest.TestLoader().loadTestsFromTestCase(TestSequencerRegress)
    unittest.TextTestRunner().run(suite)
