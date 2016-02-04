'''
Created on Jan 8, 2016

@author: dulupinar
'''
import matplotlib.pyplot as plt
import unittest
import difflib
import sequencingtools.tools as tools
import sys
import glob
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
        refReadFile = "/Users/dulupinar/Documents/UCLA/Classes/Winter16/CM222/HW2/practice_W_1/ref_practice_W_1_chr_1.txt"
        readPrefix = "/Users/dulupinar/Documents/UCLA/Classes/Winter16/CM222/HW2/practice_W_1/reads/splitreads*"
        gen_file = "../HW1/forcredit/testcase_gen.txt"
        solution_file = "../HW1/forcredit/testcase.txt"
        output = open(gen_file,'w') 
        
        readFiles = [file for file in glob.glob(readPrefix)]
        
        refSeq = tools.ReferenceSequence(refReadFile)
        kmerMap = tools.generateKmerMap(refSeq.refRead,50)
        read_count = 0
        for readFile in readFiles:
            unmatchedReads = tools.readRead(readFile)    
            len_unmatchedReads = len(unmatchedReads)
           
            i = 0
        
            for read in unmatchedReads:
                #read is a tuple
                i+=1
                
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
                        
                    refSeq.findInDels(check_read,start + 50,start + 250) #need to check 50 after start
            read_count+=1
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
