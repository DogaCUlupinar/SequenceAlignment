'''
Created on Feb 13, 2016

@author: dulupinar
'''
from collections import defaultdict, Counter
from tools import readRead
import numpy as np
import matplotlib.pyplot as plt

MIN_SIZE = 35 #try 37 too
def my_de_bruijn(sequence_reads, k):
    #sequence reads is a tuple and makes it into k
    
    de_bruijn_counter = defaultdict(Counter)
    # You may also want to check the in-degree and out-degree of each node
    # to help you find the beginnning and end of the sequence.
    for read_tuple in sequence_reads:
        
        #read the tuple into k-d mers
        pair_one = read_tuple[0]
        pair_two = read_tuple[1]
        
        for i in range(len(pair_one) -k):
            first_pattern  = pair_one[i:i+k]
            second_pattern = pair_two[i:i+k]
            node_prefix = first_pattern[:-1] + "|"+ second_pattern[:-1]
            node_suffix = first_pattern[1::] + "|"+ second_pattern[1::] 
            de_bruijn_counter[node_prefix].update([node_suffix])

    # This line removes the nodes from the DeBruijn Graph that we have not seen enough.
    de_bruijn_graph = {key: {val for val in de_bruijn_counter[key] if de_bruijn_counter[key][val] > 1}
                       for key in de_bruijn_counter}

    # This line removes the empty nodes from the DeBruijn graph
    de_bruijn_graph = {key: de_bruijn_graph[key] for key in de_bruijn_graph if de_bruijn_graph[key]}
    return de_bruijn_graph

def simple_de_bruijn(sequence_reads, k):
    """
    Creates A simple DeBruijn Graph with nodes that correspond to k-mers of size k.
    :param sequence_reads: A list of reads from the genome
    :param k: The length of the k-mers that are used as nodes of the DeBruijn graph
    :return: A DeBruijn graph where the keys are k-mers and the values are the set
                of k-mers that
    """
    de_bruijn_counter = defaultdict(Counter)
    de_bruijn_degree  = Counter()
    # You may also want to check the in-degree and out-degree of each node
    # to help you find the beginnning and end of the sequence.
    for paired_read in sequence_reads:
        for read in paired_read:
            # Cut the read into k-mers
            kmers = [read[i: i + k] for i in range(len(read) - k)]
            for i in range(len(kmers) - 1):
                pvs_kmer = kmers[i]
                next_kmer = kmers[i + 1]
                de_bruijn_counter[pvs_kmer].update([next_kmer])
                

    # This line removes the nodes from the DeBruijn Graph that we have not seen enough.
    de_bruijn_graph = {key: {val for val in de_bruijn_counter[key] if de_bruijn_counter[key][val] > 1}
                       for key in de_bruijn_counter}

    # This line removes the empty nodes from the DeBruijn graph
    de_bruijn_graph = {key: de_bruijn_graph[key] for key in de_bruijn_graph if de_bruijn_graph[key]}
    
    return de_bruijn_graph

def allDegrees(de_bruijn_graph):
    degree_dict = dict()
    for k in de_bruijn_graph:
        in_degree = inDegree(k, de_bruijn_graph)
        out_degree = len(de_bruijn_graph[k])
        degree_dict[k] = (in_degree,out_degree)
        
    return degree_dict

def de_bruijn_reassemble_smart(de_bruijn_graph,degree_dict):
    """
    Traverses the DeBruijn Graph created by simple_de_bruijn and
    returns contigs that come from it.
    :param de_bruijn_graph: A De Bruijn Graph
    :return: a list of the
    """

    assembled_strings = []
    i = 0
    skip = False
    while True:
        n_values = sum([len(de_bruijn_graph[k]) for k in de_bruijn_graph])
        print "values left ",n_values
        if n_values == 0:
            break
        found = False
        #first check for 
        for key_node in degree_dict:
            out_degree = degree_dict[key_node][1]
            in_degree  = degree_dict[key_node][0]
            if in_degree == 0 and out_degree != 0:
                good_start = key_node
                found = True
                print "found best",good_start,i                
                break
            
        if found == False:
            for key_node in degree_dict:
                out_degree = degree_dict[key_node][1]
                in_degree  = degree_dict[key_node][0]
                if in_degree < out_degree:
                    good_start = key_node
                    found = True
                    print "found better",good_start,i                
                    break
            
        if found == False:
            print "DID NOT FIND",i
            good_starts = [k for k in de_bruijn_graph if de_bruijn_graph[k]]
            good_start = good_starts[0]
            print "DID NOT FIND using",i,good_start
        # You may want to find a better start
        # position by looking at in and out-degrees,
        # but this will work okay.
        i+=1
        current_point = good_start
        assembled_string = current_point
        contig_size = 0
        while True:
            try:
                next_values = de_bruijn_graph[current_point]
                next_edge = None
                for value in next_values:
                    #find nodes where indegree = outdegree
                    try:
                        if degree_dict[value][0] == degree_dict[value][0]:
                            next_edge = value
                            next_values.remove(value)
                            print "Found good next node"
                            break
                    except KeyError:
                        break
                if next_edge ==None : next_values.pop()
                #update incoming
                degree_dict[next_edge] = (max(0,degree_dict[next_edge][0] -1),degree_dict[next_edge][1])
                #update outgoing
                degree_dict[current_point] = (degree_dict[current_point][0],max(0,degree_dict[current_point][1] -1))
                assembled_string += next_edge[-1]
                de_bruijn_graph[current_point] = next_values
                current_point = next_edge
                contig_size+=1
            except KeyError:
                print "contig Size; ",contig_size
                
                print "key error,",current_point,max(0,degree_dict[current_point][1] -1)
                degree_dict[current_point] = (degree_dict[current_point][0],max(0,degree_dict[current_point][1] -1))
                if contig_size > MIN_SIZE: assembled_strings.append(assembled_string)
                contig_size = 0
                break
    return assembled_strings

def inDegree(kmer,de_bruijn_graph):
    in_degree = 0
    for k in de_bruijn_graph:
        if kmer in de_bruijn_graph[k]:
            in_degree+=1
    return in_degree 

def createDistribution(de_bruijn_graph):
    coverage = Counter()
    
    """
    for paired_end_read in paired_end_reads:
        for read in paired_end_read:
            coverage[read]+=1
    """
    for key in de_bruijn_graph:
        for edges in de_bruijn_graph[key]:
            coverage[edges]+=1
    
    max_occurence = coverage.most_common(1)[0][1]
    distrib = np.zeros(max_occurence + 1,dtype=int)
    for key in coverage:
        distrib[coverage[key]]+=1
        
    return distrib
    
if __name__ == "__main__":
    reads = "/Users/dulupinar/Documents/UCLA/Classes/Winter16/CM222/HW3/spectrum_A_1/reads_spectrum_A_1_chr_1.txt"
    reads = "/Users/dulupinar/Documents/UCLA/Classes/Winter16/CM222/HW3/reads_hw3all_A_3_chr_1.txt"
    
    paired_end_reads = readRead(reads)
    simple_graph= simple_de_bruijn(paired_end_reads,35)
    distrib = createDistribution(simple_graph)
    plt.plot(distrib)
    plt.show()
if False:

    print "Created Graph"
    
    degree_dict = allDegrees(simple_graph)
    print "Created Dict"
    output = de_bruijn_reassemble_smart(simple_graph,degree_dict)
    #output = de_bruijn_reassemble(simple_graph)
    chr_name = "hw3all_A_3_chr_1"
    output_fn = '/Users/dulupinar/Documents/UCLA/Classes/Winter16/CM222/HW3/assembled_{0}_Super_smart_{1}.txt'.format(str(MIN_SIZE),chr_name)
    
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + chr_name + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(output))
    
