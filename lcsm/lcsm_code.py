#problem 5 @comp383
#Finding a Shared Motif

from Bio import SeqIO                                         #reads in fasta files
import os

mydir =  os.getcwd()                                           #directory
myhandle = open(mydir + "/lcsm_dataset.txt")                   #full infile directory


def commonSub(handle):
    strings = []
    records = list(SeqIO.parse(handle, "fasta"))              #make a list of SeqRecord enteries
    for record in records:                                    #from each record obtain sequence and append to strings list
        strings.append(record.seq)
    strlen = len(strings)                                     #length of strings list
    s = strings[-1]                                           #assign last string as s
    slen = len(s)                                             #length of string s

    longest = ""                                              #initialize longest common substring

    for i in range(slen):
        for j in range(i + 1, slen + 1):
            sub = s[i:j]
            k = 1
            for k in range(1, strlen):                        #for all strings in strings list
                if sub not in strings[k]:                     #if substring is not in any given string
                    break                                     #break, compare next substring

            if(k + 1 == strlen and len(longest) < len(sub)):  #if substring is found in all strings and is longer than current longest substring
                longest = sub                                 #assign substring as longest

    return longest

print(commonSub(myhandle))