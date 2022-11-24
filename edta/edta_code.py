#problem 10 @ comp 383
#Edit Distance Alignment

from Bio import SeqIO                                                   #reads in fasta files
import os


mydir =  os.getcwd()                                          #directory
myhandle = open(mydir + "/edta_dataset.txt")                  #full infile directory




def traceback(mat, s, t):
    S = ''                                                              #initialize empty string S (s')
    T = ''                                                              #initialize empty string T (t')
    n = len(s)                                                          #length of s
    m = len(t)                                                          #length of t
    i = n                                                               #intial i index for matrix (start at lowest + rightest cell)
    j = m                                                               #intial j index for matrix
    inds = n - 1                                                        #initial index for s
    indt = m - 1                                                        #initial index for t
    while(mat[i][j] != ''):                                             #while cell not empty (not top left cell)

        if(mat[i][j] == "d"):                                           #if min score obtained from diagnonal
            i = i - 1                                                   #move to prev diagonal cell (diagonal)
            j = j - 1
            S = s[inds] + S                                              #build alignment from right to left
            T = t[indt] + T
            inds = inds - 1                                              #update s, t index
            indt = indt - 1

        elif(mat[i][j] == "h"):                                          #if min score obtained from horizontal
            j = j - 1                                                    #move to prev horizontal cell (left)
            S = "-" + S                                                  #add gap "-" to sequence S
            T = t[indt] + T
            indt = indt - 1

        elif(mat[i][j] == "v"):                                          #if min score obtained from vertical
            i = i - 1                                                    #move to prev vertical cell (up)
            T = "-" + T                                                  #add gap "-" to sequence T
            S = s[inds] + S
            inds = inds - 1

    return S +  '\n'+ \
           T


def align(handle):
    records = list(SeqIO.parse(handle, "fasta"))                        #list of SeqRecord entries
    s = records[0].seq                                                  #set s as protein string taken from list of SeqRecord enteries
    t = records[1].seq                                                  #set t as protein string taken from list of SeqRecord enteries
    n = len(s)                                                          #length of s
    m = len(t)                                                          #length of t
    mat = [[0 for i in range(m + 1)] for j in range(n + 1)]             #initiliaze matrix of size (n+1)*(m+1) with 0s used to keep track of edit distance
    back = [["" for i in range(m + 1)] for j in range(n + 1)]           #initialize matrix of size (n+1)*(m+1) withs ''s used to backtrack

    #initialize ed matrix top row and first column 1...n+1, 1...m+1
    for i in range(1, n + 1):
        mat[i][0] = i
    for j in range(1, m + 1):
        mat[0][j] = j

    #initialize backtrack matrix top row and first column 'h', 'v'
    for i in range(1, n + 1):
        back[i][0] = "v"
    for j in range(1, m + 1):
        back[0][j] = "h"

    #calculate edit distance + fill in matrix
    for i in range(1, n + 1):                                           #for each cell in matrix
        for j in range(1, m + 1):

            if s[i - 1] == t[j - 1]:                                    #determine if match or mismatch, ms, by comparing
                ms = 0

            else:
                ms = 1

            p = 1                                                       #gap, p
            diagonal = mat[i - 1][j - 1] + ms                           #diagonal score
            horizontal = mat[i][j - 1] + p                              #horizontal score
            vertical = mat[i - 1][j] + p                                #vertical score

            mat[i][j] = min(diagonal,                                   #calculate score for each operation deletion, insertion, replacement) by iterativley adding previously calculated values in matrix
                            horizontal,                                 #take minimum value from each score calculation vertical, horizontal, and diagonal
                            vertical)

            #determine 'd', 'h', or 'v' + fill in matrix
            if(min(diagonal, horizontal, vertical) == diagonal):
                back[i][j] = "d"
            elif (min(diagonal, horizontal, vertical) == horizontal):
                back[i][j] = "h"
            elif (min(diagonal, horizontal, vertical) == vertical):
                back[i][j] = "v"


    distance = str(mat[n][m])                                          #lowest + rightest cell is edit distance
    alignment = traceback(back, s, t)                                  #call traceback function

    return distance + '\n' + \
           alignment                                                   #format edit distance and alignment




print(align(myhandle))