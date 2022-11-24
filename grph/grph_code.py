#problem 12 @ comp 383
#Overlap Graphs

from Bio import SeqIO
import itertools
import os

mydir =  os.getcwd()                                          #directory
myhandle = open(mydir + "/grph_dataset.txt")                   #full infile directory


def overlap(handle):
    k = 3                                                      #positive int k = 3
    datadict = {}                                              #empty dict to store record id and record seq
    records = list(SeqIO.parse(handle, "fasta"))               #list of SeqRecord entries
    for record in records:
        datadict[record.id] = record.seq                       #add record id keys, record seq values
    print(datadict)
    O = []                                                     #empty adjacency list
    for k1, k2 in itertools.combinations(datadict, 2):         #combinations method from itertools module pulls all pairwise (2) keys from dict
        s = datadict[k1]                                       #assign key 1 as s
        t = datadict[k2]                                       #assign key 2 as t
        if(s[-k:] == t[:k]):                                   #if overlap present for k1, k2 where suffix of s = prefix of t
            O.append(k1 + " " + k2)                            #append to adjacency list
        if(t[-k:] == s[:k]):                                   #if overlap present for k1, k2 where suffix of t = prefix of s
            O.append(k2 + " " + k1)                            #append to adjacency list

    for overlap in O:
        print(overlap)

overlap(myhandle)