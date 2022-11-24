#problem 9 @ comp 383
#Edit Distance

from Bio import SeqIO                                                   #reads in fasta files
import os


mydir =  os.getcwd()                                          #directory
myhandle = open(mydir + "/edit_dataset.txt")                  #full infile directory


def editdist(handle):
    records = list(SeqIO.parse(handle, "fasta"))                        #list of SeqRecord entries
    s = records[0].seq                                                  #set s as protein string taken from list of SeqRecord enteries
    t = records[1].seq                                                  #set t as protein string taken from list of SeqRecord enteries
    n = len(s)                                                          #length of s
    m = len(t)                                                          #length of t

    mat = [[0 for i in range(m + 1)] for j in range(n + 1)]            #initiliaze matrix of size (n+1)*(m+1)

                                                                       #initialize top row and first column 1...n+1, 1...m+1
    for i in range(1, n + 1):
        mat[i][0] = i
    for j in range(1, m + 1):
        mat[0][j] = j

    for i in range(1, n + 1):                                          #for each cell in matrix, starting with top left cell
        for j in range(1, m + 1):

            if s[i - 1] == t[j - 1]:                                   #determine if match or mismatch, ms, by comparing
                ms = 0

            else:
                ms = 1

            p = 1                                                      #gap, p

            mat[i][j] = min(mat[i - 1][j] + p,                         #calculate score for each operation deletion, insertion, substitution by iterativley adding previously calculated values in matrix
                            mat[i][j - 1] + p,                         #take minimum value from each score calculation vertical, horizontal, and diagonal
                            mat[i - 1][j - 1] + ms)

    distance = mat[n][m]                                               #lowest + rightest cell is edit distance
    return distance




print(editdist(myhandle))