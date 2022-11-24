#problem 8 @ comp 383
#Introduction to Set Operations

import os

mydir =  os.getcwd()                                           #directory
myhandle = open(mydir + "/seto_dataset.txt")                   #full infile directory

def sets(handle):
    items = handle.readlines()                                              #list of lined read from files
    sets = []                                                               #empty dictionary sets
    for item in items:                                                      #for each line read from file
        sets.append(item.rstrip("\n"))                                      #append to sets, remove "\n"
    n = int(sets[0])                                                        #set n as positive integer
    U = {i for i in range(1,n + 1)}                                         #create set U consisting of {1,2,…,n}
    a = sets[1].replace("{", "").replace("}", "").replace(" ", "")          #reformat subset str, remove '{', '}', and '
    b = sets[2].replace("{", "").replace("}", "").replace(" ", "")
    lista, listb = a.split(","), b.split(",")                               #create a list of set values by splitting subset strs at commas
    lista, listb = [int(i) for i in lista],[int(i) for i in listb]          #convert list contents from str to int
    A, B = set(lista), set(listb)                                           #convert subset lists to subsets
    union = A | B                                                           #union
    intersection = A & B                                                    #intersection
    diffA = A - B                                                           #set difference
    diffB = B - A                                                           #set difference
    setcompA = U - A                                                        #set complement, Ac
    setcompB = U - B                                                        #set complement, Bc
    print(union)
    print(intersection)
    print(diffA)
    print(diffB)
    print(setcompA)
    print(setcompB)

sets(myhandle)