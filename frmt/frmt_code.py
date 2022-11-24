#problem 2 @ comp 383
#Data Formats


from Bio import Entrez
from Bio import SeqIO
import os

mydir =  os.getcwd()                                          #directory
handle = open(mydir + "/frmt_dataset.txt")                    #handle object
ids = handle.read()                                           #GenBank entry objects

Entrez.email = "apatel72@luc.edu"


def shortest_seq(IDs):
    handle = Entrez.efetch(db="nucleotide", id=[IDs], rettype="fasta")  #obtains plain text records in FASTA format from NCBI's [Nucleotide] database
    records = list(SeqIO.parse(handle, "fasta"))            #list of SeqRecord entries
    shortseq = float("inf")                                 #intialize variable
    for record in records:                                  #interate through records
        if shortseq > len(record.seq):                      #if length of current record seq is smaller set as shortseq
            shortseq = len(record.seq)
            shortest = record
    print(shortest.format("fasta"))                         #shortest record in FASTA format


shortest_seq(ids)