#problem 15 @ comp 383
#Genome Assembly Using Reads

from Bio.Seq import Seq
import itertools
import sys
import os
sys.setrecursionlimit(5000)                                                                        #recursion limit was needed to be reset because of recursive detect_cycle/size of data set
                                                                                                   #better method would be to rewrite recursive function as an iterative function
                                                                                                   #recursion limit does not need to be reset for small sample data set

mydir =  os.getcwd()                                          #directory
myhandle = open(mydir + "/gasm_dataset.txt")                  #full infile directory


temp = []                                                                                          #list to store remaining nodes to be visited
visited = set()                                                                                    #set of adjacent vertices to keep track
assemblies = []                                                                                    #list of super created per cycle 

def detect_cycle(visited,graph,node,i, assembly):                                                  #recursive func to detect cycle 
    if (node not in graph.keys()):                                                                 #if node not found as key in graph determine no Cycle
        return "No Cycle", assemblies
    elif (node not in visited):                                                                    #if node not visited
        visited.add(node)                                                                          #visit node + add to visited set

        for neighbor in graph[node]:                                                               #for an adjacent neighbor node whose prefix overlaps current node
            assembly = assembly + neighbor[-1]                                                     #append last letter of adj to assembly
            return detect_cycle(visited, graph, neighbor, i, assembly)                             #call function again with neighbor node as new node to be explored

    elif (node in visited):                                                                        #if node has already been explored
        assemblies.append(assembly)                                                                #append superstring of cycle to assemblies
        i = i + 1                                                                                  #add 1 to i to keep track of cycle count
        temp = [graph[node] for node in graph if node not in visited]                              #create a list of nodes that are in graph but still need to be visited
        flat_temp = list(itertools.chain(*temp))

        if (len(temp) == 0):                                                                       #if there are no more such nodes in temp
            return str(i) + " Cycle(s)", assemblies[0]                                             #return number of cycles, and one superstring


        else:                                                                                      #if there are still nodes that need to be visited
            node = flat_temp[0]                                                                    #pick a node from temp, set as new node
            assembly = ""                                                                          #restart superstring for new cycle
            for neighbor in graph[node]:
                assembly = assembly + neighbor[-1]                                                 #continue building superstring and calling function
                return detect_cycle(visited, graph, neighbor,i,assembly)




def what(handle):
    items = handle.readlines()                                                                    #list of lines read from files
    k1mers = []                                                                                   #empty dictionary for k1mers
    S = set()                                                                                     #set for (k+1)-mers
    Src = set()                                                                                   #set for reverse complement(k+1)-mers
    for item in items:                                                                            #for each line read from file
        k1mers.append(Seq(item.rstrip("\n")))                                                     #append each (k+1)-mer to k1mers list, remove "\n"
    [S.add(x) for x in k1mers]                                                                    #create a set S of non duplicate (k+1)-mers
    [Src.add(x.reverse_complement()) for x in S]                                                  #create a set Src of non duplicate reverse complement (k+1)-mers
    u = S | Src                                                                                   #S∪Src, union
    temp_u = set()                                                                                #temporary set u
    cont = True

    while(cont != False):                                                                         #while continue is not False or until 2 cycles are detected
        edges = []                                                                                #list of edges that will update every loop
        G = {}                                                                                    #empty graph dict
        temp_u.clear()                                                                            #clear temp_u before starting next loop
        for k1mer in u:                                                                           #for (k+1)-mer in set U
            k = len(k1mer) - 1                                                                    #set value of k (k = len(k+1)-mer) - 1
            s = str(k1mer[0:k])                                                                   #(k+1)-mer's kmers present as substring
            t = str(k1mer[1:k + 1])
            temp_u.add(s)                                                                         #add new kmers to temp_u
            temp_u.add(t)

            G[s] = [t]                                                                            #create graph
            edges.append((s, t))                                                                  #append kmers edges as tuples to adj list

        edges = list(set(edges))                                                                  #remove duplicates in case there are any

        u.clear()                                                                                 #clear u
        [u.add(k1mer) for k1mer in temp_u]                                                        #update u with new kmer values from temp_u

        node = edges[0][0].replace("'", "")                                                       #pick a kmer from edges as a node

        search = detect_cycle(visited,G,node,0,'')                                                #check for cycle, call function detect_cycle
        cycle_detected = search[0]                                                                #first return value
        superstring = search[1]                                                                   #second return value, superstring that was built in detect_cycle function

        if(cycle_detected == "2 Cycle(s)"):                                                       #if 2 cycles are detected
            cont = False                                                                          #set cont to False and exit while loop
        else:
            cont = True                                                                           #else continue

    return superstring                                                                            #return superstring




print(what(myhandle))