#problem 3 @comp383
#Transitions and Transversions


from Bio import SeqIO                                         #reads in fasta files
import os

mydir =  os.getcwd()                                           #directory
myhandle = open(mydir + "/tran_dataset.txt")                  #full infile directory

def tt(handle):
    records = list(SeqIO.parse(handle, "fasta"))              #make a list of SeqRecord enteries
    seq1 = records[0].seq                                     #obtain string 1 and string 2 from records
    seq2 = records[1].seq
    seqlen = len(seq1)                                        #string length
    transversions = 0
    transitions = 0
    for x in range(0, seqlen-1):                              #for every position in both strings compare bases, for a transition or transversion add 1
        if(seq1[x] == seq2[x]):
            transitions = transitions + 0
            transversions = transversions + 0
        elif(seq1[x] == "A" and seq2[x] == "G"):
            transitions = transitions + 1
        elif(seq1[x] == "A" and seq2[x] == "C"):
            transversions = transversions + 1
        elif(seq1[x] == "A" and seq2[x] == "T"):
            transversions = transversions + 1
        elif(seq1[x] == "G" and seq2[x] == "A"):
            transitions = transitions + 1
        elif(seq1[x] == "G" and seq2[x] == "C"):
            transversions = transversions + 1
        elif(seq1[x] == "G" and seq2[x] == "T"):
            transversions = transversions + 1
        elif(seq1[x] == "C" and seq2[x] == "T"):
            transitions = transitions + 1
        elif(seq1[x] == "C" and seq2[x] == "A"):
            transversions = transversions + 1
        elif (seq1[x] == "C" and seq2[x] == "G"):
            transversions = transversions + 1
        elif(seq1[x] == "T" and seq2[x] == "C"):
            transitions = transitions + 1
        elif(seq1[x] == "T" and seq2[x] == "A"):
            transversions = transversions + 1
        elif (seq1[x] == "T" and seq2[x] == "G"):
            transversions = transversions + 1

    ratio = transitions/transversions                          #calculate ratio
    return ratio

print(tt(myhandle))