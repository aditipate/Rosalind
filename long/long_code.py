#problem 14 @ comp 383
#Genome Assembly as Shortest Superstring

from Bio import SeqIO
import itertools
import os

mydir =  os.getcwd()                                           #directory
myhandle = open(mydir + "/long_dataset.txt")                   #full infile directory


def update_edges(edges,s,t):
    edges = list(filter(lambda x: x[0] != s, edges))                                               #update edges, filter out edges containing vertices that have been merged
    edges = list(filter(lambda x: x[1] != t, edges))
    edges = list(filter(lambda x: x[0] != t, edges))
    edges = list(filter(lambda x: x[1] != s, edges))

    return edges

def update_vertices(vertices,s,t,o):
    merge = s + t[o:len(t) + 1]                                                                    #create new merged vertex
    vertices.remove(s)                                                                             #remove vertices that were merged
    vertices.remove(t)
    vertices.append(merge)                                                                         #append new merged vertex to vertices

    return vertices


def scs(handle):
    vertices = []                                                                                   #empty list to store vertices
    edges = []                                                                                      #empty list to store edges

    records = list(SeqIO.parse(handle, "fasta"))                                                    #list of SeqRecord entries
    for record in records:
        vertices.append(record.seq)                                                                 #store record seq enteries as vertices

    while(len(vertices) > 1):                                                                       #while there are still vertices to be merged
        for s, t in itertools.combinations(vertices, 2):                                            #combinations method from itertools module pulls all pairwise (2) keys from dict
                curr = 0                                                                            #curr is the index that currently gives the longest/maximal overlap bwtween 2 keys s,t
                order = ""                                                                          #order determines which key, s or t, generates suffix or prefix
                for k in range(1, len(s)):                                                          #for all possible kmer lengths k in s and t
                    if (s[-k:] == t[:k]):                                                           #check to see if overlap present between suffix of s and prefix of t
                        if (k > curr):                                                              #if overlap given by k is longer than current length of over lap, set curr = k
                            curr = k
                            order = "s,t"                                                           #set appropriate order
                    elif (t[-k:] == s[:k]):                                                         #check to see if overlap present between suffix of t and prefix of s
                        if (k > curr):                                                              #if overlap given by k is longer than current length of over lap, set curr = k
                            curr = k
                            order = "t,s"                                                           #set appropriate order

                if (order == "s,t"):                                                                #based on order append s,t pair and curr to edges
                    edges.append((s, t, curr))
                elif (order == "t,s"):
                    edges.append((t, s, curr))

        longest = max(edges, key = lambda item: item[2])                                            #find the edge with the greatest curr value or longest overlap associated with it
        s = longest[0]                                                                              #retrieve reads to be merged
        t = longest[1]
        o = longest[2]
        edges = update_edges(edges,s,t)                                                             #update edges
        vertices = update_vertices(vertices,s,t,o)                                                  #update vertices, adding newly merged vertex

    shortest_super = vertices[0]                                                                    #return last remaining vertex in vertices or shortest superstring
    print(shortest_super)


scs(myhandle)