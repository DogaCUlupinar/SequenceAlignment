from collections import defaultdict
from distutils.sysconfig import PREFIX

def CyclicChromosomes(reads):
    #reads is a list of reads
    cyclic_string = ""
    read = reads.pop()
    while reads:
        suffix = read[1:len(read)]
        for pre in reads:
            prefix = pre[:len(pre)-1]
            if suffix == prefix:
                cyclic_string+=suffix[-1]
                break
        read = pre
        reads.remove(pre)
    cyclic_string+=pre[-1]
    return cyclic_string

def readReads(filename):
    reads = []
    with open(filename) as f:
        for line in f:
            reads.append(line.strip())
            
    return reads

def debruijnGraph(reads):
    #reads is a list of reads
    graph = defaultdict(list)
    for read in reads:
        prefix = read[0:len(read)-1]
        suffix = read[1:len(read)]
        graph[prefix].append(suffix)
    
    return graph

def createAL(graph):   
    #create the adjaceny list
    adjaceny_list = ""
    for node in graph:
        adjaceny_list+="{suff} -> {edges}\n".format(suff=node,edges=",".join(graph[node])) 
    
    return adjaceny_list

def inDegree(kmer,de_bruijn_graph):
    in_degree = 0
    for k in de_bruijn_graph:
        in_degree+=de_bruijn_graph[k].count(kmer)
    return in_degree 

def allDegrees(de_bruijn_graph):
    degree_dict = dict()
    for k in de_bruijn_graph:
        in_degree = inDegree(k, de_bruijn_graph)
        out_degree = len(de_bruijn_graph[k])
        degree_dict[k] = (in_degree,out_degree)
        
        for edge in de_bruijn_graph[k]:
            if edge not in de_bruijn_graph:
                in_degree = inDegree(edge, de_bruijn_graph)
                degree_dict[edge] = (in_degree,0)
        
    return degree_dict

def runProblem4():
    reads = ["ATG","ATG","TGT","TGG","CAT","GGA","GAT","AGA"]
    reads = readReads("rosalind_ba3k.txt")
    graph = debruijnGraph(reads)
    degree_dict = allDegrees(graph)
    first_degree = allDegrees(graph)
    contigs = []
    
    while sum([len(graph[node]) for node in graph]) >0:
        good_starts = [key for key in degree_dict if degree_dict[key] != (1,1) and degree_dict[key][1] > 0]
        good_start = good_starts[0] if len(good_starts) else degree_dict.itervalues().next()[0]
        
        contig = good_start
        dont_stop = True
        while(dont_stop):
            dont_stop = False
            possible_edges = graph[good_start]
            if len(possible_edges) == 0: break
            for edge in possible_edges:
                if first_degree[edge] == (1,1):
                    dont_stop = True
            contig+=edge[-1]
            
            #update incoming
            degree_dict[edge] = (max(0,degree_dict[edge][0] -1),degree_dict[edge][1])
            #update outgoing
            degree_dict[good_start] = (degree_dict[good_start][0],max(0,degree_dict[good_start][1] -1))
                
            possible_edges.remove(edge)
            good_start = edge
        
        contigs.append(contig)
    return "\n".join(contigs)
    
def runProblem1():
    test_reads = ["ATTAC","TACAG","GATTA","ACAGA","CAGAT","TTACA","AGATT"]
    filename = "rosalind_pcov.txt"
    reads = readReads(filename)
    print CyclicChromosomes(reads)
    
def runProblem2():
    test_reads = ["GAGG","CAGG","GGGG","GGGA","CAGG","AGGG","GGAG"]
    reads = readReads("rosalind_ba3e.txt")
    print createAL(debruijnGraph(reads))

def runProblem3():
    dna = "AGCCATTAGAAGAGACGTGGTACATGGTCCCAATCATCCACGGGCTTTTGGTGCGAGCTCGTAATAATCTCATACTTCTACAACCCCCTCACGATTCTAAGCGTTGCAAATTGGAGTAACGACGTCCGATAGCAGTTATCCGTACCGCCGTAGATTTCAATCACGTTGACTCAACAATTTGTTCCTCTACAGTCTCTTAGAAGAGGACAGACCAAGGTAATACGGCGCGCCGCAAAAGTCGCCTGCTCGAACACCCGTAAAGAGTGAGGGATCATGTAGTAGACGCCTAAACACAACGGCTGTCACTGTGCAGCTGCCGAATGTTGACCACAGTCTGTAAGGAACCGAAGGGTGTCCGTGGCGTTGGACTGAGCAGTAACCCTGTTCGTTTACACTTTGAGGAAAGACCGTCTGCTGTAGAGCTCCTGCTCTGAAATTGAGTTGCACAAAAAATGATTGTGACTGTATGCGAACTGTGACTTCATGAAGAATCTCCTAAACATCACAAGATTGAGACGCGAGCCCAAAACCCTGCTTCGACACGATGCTACGGCGACATTGAGAACATTCGCCCTATGTTGCTGAGGTGTGGACACACAGCCACGCGGTGTAAATCTATGGTCTCGAATAGGAGTGGCAGCGAATAATACGGGCTGTATTATCTAAAGAAATGAGGAGCACTCCCGAGCTACAATCATGCACATTGACGGACAATACATCCTACGGGTGGAACTTGCCCTATGTAAGGAACCGATTCAATCTACCAATTGTGTACCGATAAGTATCGCCCCGCCCTAGACTAGGAGGGGCGTGCCTCCGTTTATGGGTGGGCCACAGATGCTTGTGCCACGCTCGGTTGCACGACTCCCTTCTGTCAGAAAGGTAATGTAGCCAACTCCGAAGGATTTGACCATACTCTGTCAAGCTCCCGATTGGACATTGTAGTATGATCAGGACCTCATTGTTATTCTCGTTAGTGTAATCAGAAATACTTGGCAGAAAACATTCTAATACTTTTAGGAAGCATGTAACACACCACCTCTCTTGATAAGAGGAGACGACAACTCGGAGGCACCAATTACGAGGTGGAGCCAGTCCGACGTGGACTTCAAGCGCTACATACTTCGATTCCGCATGGTTTTTCAATAGTACCAACGGATCATTATTAAGTAAATACTTTAAATAGGCTCTCAGAACACGGTAGGCGCGCATGTCGTATCGGTCCACCCACCCCGGCAAGTCGGGATAGATGGTAAGGCCATGCCACGCATGCAGTAAATTCTCTATCTTAACTCTAAGGGCTTCATCCTGTCTTAGACGGGCCCATGGAGGGGCAGGTTTTGTCGATTGAGAACAGTGCGCCGGGTTTTTGCTCCGTATCGGTCTAGATTATTTATGCGTTGAGACTTCCAATCGCAGTGCTCTCTGGGAGTAACAGGGATTCGCCTAAGCCGGGTATAGCATGTTTCTAACCGCATAAAACTTTAGTCATGCAGTGACATTAGGCAACTCGAATTATAGCCTTGATGTCACTAACCGGAAAGGTGGCTGGAGCCCCGAAGCCTTCTGCTTAATGTTTACTCAAGCTCACCGGACATCCACCAGATACCAGTGGTGGCCCAAGACGACGGGATGTAGCTCGAGCCTTAATATTAGTGGGGAATTGCGCCCACGCCATCAACTAGCAACAAGCTATGACAGGTCCCTCGTTCTGTACTACCTGGTGCTATTCCTTTTGCAAAATTTCTGAAGGTGCAGCTAGCTCAATGGAGGCCCCGTAACACTTACTATTGTACCTACATGCTAGATACTACTACCCAATCCTTTTCCGGGAACATCGAGCGGACCATTGGCGGCTATCCCTAAATGACTTTATAGATTACGGTGGATCTAATGTGACCCTTTTGTTGTTTAAAGTAGACAATTGGCATAATAAATTGATACCCCCACACCCGGCAATGCTCCATTCTGAGATATATTAGGCCGCCAGTGTATGTTTC"
    kmer_len = 12
    reads = []
    for i in range(len(dna)-kmer_len +1):
        reads.append(dna[i:i+kmer_len])
    
    print createAL(debruijnGraph(reads))
   
if __name__ == "__main__":
    print "START"
    print runProblem4()

    