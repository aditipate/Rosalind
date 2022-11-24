#problem 6 @ comp 383
#RNA Splicing

from Bio.Seq import Seq                                       #Seq methods
from Bio import SeqIO                                         #reads in fasta files
import os

mydir =  os.getcwd()                                           #directory
myhandle = open(mydir + "/splc_dataset.txt")                   #full infile directory


def rna_splice(handle):
    records = list(SeqIO.parse(handle, "fasta"))              #list of SeqRecord entries
    substrs = []                                              #empty list for substrings
    dnastr = str(records[0].seq)                              #dna string taken from list SeqRecord enteries, convert to str object in order to use replace method
    for i in range(1,len(records)):                           #list of remaning strings or substrings acting as introns from SeqRecord enteries
        substrs.append(str(records[i].seq))
    for substr in substrs:                                    #for each substring check to see if present in dna string
        if(substr in dnastr):
            dnastr = dnastr.replace(substr, "")               #if substring present, remove instance from dna string
        else: continue
    exons = Seq(dnastr)                                       #create new Seq object from final dna string (introns removed)
    mrna = exons.transcribe()                                 #transcribe exons
    protein = mrna.translate()                                #translate mrna
    return protein[:-1]                                       #format protein sequence, remove '*'



print(rna_splice(myhandle))