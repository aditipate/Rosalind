#problem 11 @ comp 383
#GenBank Introduction

from Bio import Entrez
import os

mydir =  os.getcwd()                                          #directory
myhandle = open(mydir + "/gbk_dataset.txt")                  #full infile directory



Entrez.email = "apatel72@luc.edu"


def genus_enteries(handle):
    filedata = list(handle.readlines())                        #list of file data
    records = list(map(lambda x: x.strip(),filedata))          #strip '\n' from data in filedata
    genus = records[0]                                         #assign genus
    date_1 = records[1]                                        #assign date 1
    date_2 = records[2]                                        #assign date 2
    dates_range = date_1 + ':' + date_2
    handle = Entrez.esearch(db="nucleotide", term=genus + '[Organism]' + ' AND ' + dates_range + '[PDAT]')  #GenBank entries for given genus Nucleotide database between dates specified
    record = Entrez.read(handle)                               #store record
    print(record["Count"])                                     #return number of GenBank enteries

genus_enteries(myhandle)