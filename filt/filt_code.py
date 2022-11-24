#problem 18 @ comp 383
#Read Filtration by Quality

from Bio import SeqIO
import os


mydir =  os.getcwd()                                          #directory
infile = open(mydir + "/filt_dataset.txt")                    #full infile directory



def filter_reads(handle):
    filtered = 0                                                                         #initiatize var to add total number of filtered reads
    data = handle.readline().split()                                                     #split contents of first line into list
    q, p = int(data[0]), int(data[1])/100                                                #set threshold, percentage
    for record in SeqIO.parse(handle, "fastq"):                                          #for each record in the fastq file
        qual = record.letter_annotations["phred_quality"]                                #obtain the Phred quality scores corresponding to each base of the sequence
        below = 0                                                                        #initialize var to add the total number of scores that are below the threshold
        for score in qual:                                                               #for each score in the list of Phred quality scores
            if(score >= q):                                                              #if score is greater or equal to threshold
                below = below + 1                                                        #increment below to keep track of scores above threshold
        if((below >= len(qual)*p)):                                                      #if number of scores below threshold are atleast p% of the total scores
            filtered = filtered + 1                                                      #incremement filtered to keep track of reads


    print(filtered)



filter_reads(infile)