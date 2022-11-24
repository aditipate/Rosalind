#problem 7 @ comp 383
#k-Mer Composition

import itertools
from Bio import SeqIO                                                      #reads in fasta files
import os

mydir =  os.getcwd()                                           #directory
myhandle = open(mydir + "/kmer_dataset.txt")                   #full infile directory

def kmercomp(handle):
    records = list(SeqIO.parse(handle, "fasta"))                           #list of SeqRecord entries
    s = records[0].seq                                                     #dna string, s taken from list of SeqRecord enteries
    A = {}                                                                 #empty 4-mer dictionary
    bases = ['A', 'T', 'G', 'C']                                           #list of nucleotide bases
    kmers = [''.join(p) for p in itertools.product(bases, repeat=4)]       #list of every possible 4-mer of nucleotide bases using itertools
    sortedkmers = sorted(kmers, key=lambda x: x.lower())                   #list of sorted kmers alphabetically
    for kmer in sortedkmers:                                               #assign every 4-mer as a key in dictionary A, intialize values as 0
        A[kmer] = 0
    for i in range(0,len(s)):                                              #iterate over dna string
        if(s[i:i+4] in A):                                                 #if 4-mer substring is a key in A
            kmer = s[i:i+4]                                                #set substring as 4-mer key
            A[kmer] = s.count_overlap(kmer)                                #set key value to the number of times 4mer occurs in dna string s, include overlap count
    composition = A.values()                                               #4-mer composition
    print(*composition, sep = " ")                                         #format, remove commas and brackets

kmercomp(myhandle)