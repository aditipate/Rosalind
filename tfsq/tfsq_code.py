#problem 16 @ comp 383
#FASTQ format introduction

from Bio import SeqIO
import os

mydir =  os.getcwd()                                                                    #directory

infile = mydir + "/tfsq_dataset.txt"                                                       #full infile directory
outfile = mydir + "/rosoutput_tfsq.txt"                                                    #full outfile directory




def fastA(input_handle):                                                                 #function to open newly written fasta file
    with open(input_handle, 'r') as input:                                               #open and read infile
        print(input.read())                                                              #print file content



def fastQ_A(input_handle,output_handle):                                                 #function to convert fastq to fasta
    with open(input_handle, "r") as input:                                               #open and read infile
        with open(output_handle, "w") as output:                                         #open outfile to write
            records = SeqIO.parse(input, "fastq")                                        #parse fastq infile and place contents into records list
            SeqIO.write(records, output, "fasta")                                        #write records to outfile in fasta format

    fastA(output_handle)                                                                 #call function fastA



fastQ_A(infile,outfile)