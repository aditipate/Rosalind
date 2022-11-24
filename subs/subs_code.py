#problem 4 @comp383
#Finding a Motif in DNA

import os

mydir =  os.getcwd()                                           #directory
myhandle = open(mydir + "/subs_dataset.txt")                   #full infile directory

alist = myhandle.readlines()                                  #list of lines read from file
alist[0] = alist[0].strip()                                   #remove /n from strings
alist[1] = alist[1].strip()
s = "0" + alist[0]                                            #assign string and substring from list to variables and append space holder or "0" to shift index of string
t = "0" + alist[1]

def locations(str,sub):
    loc = []                                                  #empty list
    for i in range(1, len(str)):                              #for the length of the string, compare every substring of same length with given substring
        if(str[i:i+len(sub)-1] == sub[1:len(sub)]):
            loc.append(i)                                     #append index substring's beginning position to list
    return loc                                                #return all locations of t as a substring of  s as a list



print(locations(s,t))