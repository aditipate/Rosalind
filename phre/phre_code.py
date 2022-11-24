#problem 17 @ comp 383
#Read Quality Distribution

from Bio import SeqIO
import os

mydir =  os.getcwd()                                           #directory
infile = open(mydir + "/phre_dataset.txt")                     #full infile directory


def belowThresh(handle):
    below = 0                                                                            #initialize var to add the total number of reads reads whose average quality is below the threshold
    thres = int(handle.readline())                                                       #read in first line, set int as threshold

    for record in SeqIO.parse(handle, "fastq"):                                          #for each record in the fastq file
        qual = record.letter_annotations["phred_quality"]                                #obtain the Phred quality scores corresponding to each base of the sequence
        average_qual = sum(qual)/len(qual)                                               #add all Phred quality scores and dived by number of scores (average)
        if(average_qual < thres):                                                        #check to see if average qaulity is below threshold
            below = below + 1                                                            #if it is below thresh increment below


    print(below)




belowThresh(infile)