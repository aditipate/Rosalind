#problem 19 @ comp 383
#Base Quality Distribution

import os
from Bio import SeqIO

mydir =  os.getcwd()                                          #directory
infile = open(mydir + "/bphr_dataset.txt")                    #full infile directory


def base_qual(handle):
    below = 0                                                                           #initialize var to add the total number base quality scores
    total = 0                                                                           #intitialize var to keep track of total records
    thres = int(handle.readline())                                                      #read in first line, set int as threshold
    base_scores = []                                                                    #empty list
    for record in SeqIO.parse(handle, "fastq"):                                         #for each record in the fastq file
        total = total + 1                                                               #add total number of records
        qual = record.letter_annotations["phred_quality"]                               #obtain the Phred quality scores corresponding to each base of the sequence
        if base_scores == []:                                                           #if base scores list empty
            base_scores = [0] * len(qual)                                               #intialize list with lists containing 0 with a length corresponding to number of phred scores
        for i in range(len(qual)):                                                      #for each base score in qual
            base_scores[i] = base_scores[i] + qual[i]                                   #add score to corresponding list in base_scores list

    for i in range(0, len(base_scores)):                                                #for each base score, get the mean
        base_scores[i] = base_scores[i] / total                                         #add mean score to corresponding list in base scores
    for mean_score in base_scores:                                                      #for each mean score in base score
        if(mean_score < thres):                                                         #if mean score is less than threshold
            below = below + 1                                                           #increment below

    print(below)



base_qual(infile)