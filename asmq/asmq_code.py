#problem 21 @ comp 383
#Assessing Assembly Quality with N50 and N75

import os

mydir =  os.getcwd()                                          #directory
infile = open(mydir + "/asmq_dataset.txt")                    ##full infile directory


def nXX(contigs):
    total = sum(contigs)                                                          #get total sum of all the contig lengths or total base pairs
    n50 = 0                                                                       #set n50 to 0
    n75 = 0                                                                       #set n75 to 0
    for i in range(0,len(contigs)):                                               #for all contig lengths
        bps = sum(contigs[i:])                                                    #get the sum of contig lengths equal and rgreater to current contig length
        if(bps > (total * 0.5) and contigs[i] > n50):                             #if sum of base pairs is greater than 50% of total base pairs
            n50 = contigs[i]                                                      #set current contig length as n50
        if (bps > (total * 0.75) and contigs[i] > n75):                           #if sum of base pairs is greater than 75% of total base pairs
            n75 = contigs[i]                                                      #set current contig length as n75

    return n50, n75


def assess_assembly(handle):
    items = handle.readlines()                                                       #list of lines read from files
    contigs = []                                                                     #empty list for contigs
    contig_lens = []                                                                 #empty list for contigs lengths
    for item in items:                                                               #for each line read from file
        contigs.append((item).rstrip("\n"))                                          #remove "\n" and append contigs to contigs list
    [contig_lens.append(len(contig)) for contig in contigs]                          #append contig lengths to contig_lens list
    contig_lens.sort()                                                               #sort contigs lengths from shortest to longest

    n50 = nXX(contig_lens)[0]                                                        #call function nXX and get n50
    n75 = nXX(contig_lens)[1]                                                        #call function nXX and get n75
    print(n50, n75)



assess_assembly(infile)