#problem 1 @ comp 383
#Computing GC Content


from Bio import SeqIO                                         #reads in fasta files
import os


mydir =  os.getcwd()                                          #directory
handle = open(mydir + "/gc_dataset.txt")                       #handle object


#function computes gc content of dna sequence
def gc(sequence):
    count = sequence.count("G") + sequence.count("C")       #count total number of Gs and Cs in sequence
    content = count / len(sequence) * 100                   #calculate gc content
    return content


#function returns the highest gc content of sequences from a fasta file and the sequence id
def highestgc(filepath):
    highest = 0                                                   #intialize variables
    id = ""
    for record in SeqIO.parse(filepath, "fasta"):                 #parse each seq in fasta file
        if gc(record.seq) > highest:                              #if gc content of current record seq is greater set as highest
            highest = gc(record.seq)
            id = record.id
    return id, highest                                            #return id and highest gc


print(highestgc(handle))