#problem 13 @ comp 383
#Constructing a De Bruijn Graph

from Bio.Seq import Seq
import os

mydir =  os.getcwd()                                          #directory
myhandle = open(mydir + "/dbru_dataset.txt")                  #full infile directory



def deBruijn(handle):
    items = handle.readlines()                                   #list of lines read from files
    k1mers = []                                                  #empty dictionary for k1mers
    edges = []                                                   #empty edge adjacency list
    S = set()                                                    #set for (k+1)-mers
    Src = set()                                                  #set for reverse complement(k+1)-mers
    for item in items:                                           #for each line read from file
        k1mers.append(Seq(item.rstrip("\n")))                    #append each (k+1)-mer to k1mers list, remove "\n"
    [S.add(x) for x in k1mers]                                   #create a set S of non duplicate (k+1)-mers
    [Src.add(x.reverse_complement()) for x in S]                 #create a set Src of non duplicate reverse complement (k+1)-mers
    U = S | Src                                                  #S∪Src, union
    for k1mer in U:                                              #for (k+1)-mer in set U
        k = len(k1mer)-1                                         #set value of k (k = len(k+1)-mer) - 1
        s = str(k1mer[0:k])                                      #(k+1)-mer's kmers present as substring
        t = str(k1mer[1:k+1])
        edges.append((s, t))                                     #append kmers edges as tuples to adj list
    edges = list(set(edges))                                     #remove duplicates in case there are any
    for edge in edges:                                           #format each edge, remove "'"
        s = edge[0].replace("'", "")
        t = edge[1].replace("'", "")
        print("(" + s + ", " + t + ")")                          #format and print all edges


deBruijn(myhandle)