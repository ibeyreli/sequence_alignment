# -*- coding: utf-8 -*-
"""
CS481 Fall 2019 Homework 3
Sequence Alignment
Source Code File
Created on Thu Nov 7 18:09:11 2019
@author: İlayda Beyreli 201801130
"""
import sys
import time
from itertools import chain

def readfa2(file_name):
    # Reading .FA files into lists
    with open(file_name, "r") as fin:
        data =fin.read()#.splitlines()
    s = data.index(">")
    s = data[s+1:].index(">")
    seq1 = data[:s].splitlines()
    n1 = seq1[0]
    seq1 = list(chain.from_iterable(seq1[1:])) 
    seq2 = data[s+1:].splitlines()
    n2 = seq2[0]
    seq2 = list(chain.from_iterable(seq2[1:])) 
    return n1, seq1, n2, seq2

def writefa2(score,n1, seq1, n2, seq2,file_name):
    if len(n1) > len(n2):
        n2+=" "*(len(n1)-len(n2))
    elif len(n1) < len(n2):
        n1+=" "*(len(n2)-len(n1))
    i = 0
    j = 0
    with open(file_name, "w+") as fout:
        fout.write("Score:"+str(score)+"\n\n")
        while j < len(seq2) and i  < len(seq1):
            fout.write(n1+" "+seq1[i:i+60]+"\n")
            i+=60
            fout.write(n2+" "+seq2[j:j+60]+"\n\n")
            j+=60
    return True

def zeros(r,c):
    # Matrix as list of lists 
    return [[0 for i in range(c)] for i in range(r)]

def matprint(matrix):
    for i in range(len(matrix)):
            print(matrix[i],"\n")
    return 0

def naive_score(c1,c2,gapop):
    # Scoring Function, given in the assignment
    if c1 == c2:
        return 2 #match_premium
    elif c1 =='-' or c2 =='-':
        return gapop
    else:
        return -3 #missmatch
def argmax(a1,a2,a3,a4):
    switcher={0:'d',
            1:'u',#del
            2:'l',#in
            3:'d'}
    l = [a1,a2,a3,a4]
    m = max(a1,a2,a3,a4)
    m = l.index(m)
    return switcher.get(m)
    
def naive_alignment(seq1,seq2,gapop,local=False):
    # Needleman-Wunsh & Smith-Whitterman Naive Alignment Algorithms
    # Initialize nm table
    n = len(seq1)
    m = len(seq2)
    S = zeros(n+1,m+1)
    L = zeros(n+1,m+1)
    if local == False:
        # Initialize according to global alignment
        l = -float("inf")
        for i in range(n+1):
            S[i][0] = i*gapop
        for i in range(m+1):
            S[0][i] = i*gapop
    else:
        # Initialize according to local alignment
        l = 0
    # Fill the scoring table    
    for i in range(1,n+1):
        for j in range(1,m+1):
            match = S[i-1][j-1]+naive_score(seq1[i-1],seq2[j-1],gapop)
            del1 = S[i-1][j]+naive_score(seq1[i-1],'-',gapop)
            in1 = S[i][j-1]+naive_score('-',seq2[j-1],gapop) 
            S[i][j] = max(match,del1,in1,l)
            L[i][j] = argmax(match,del1,in1,l)
    if local == False:
        # Score for global alignment
        alignscore = S[-1][-1]
    else:
        # Score for local alignment
        alignscore = max([max(row) for row in S])
    # matprint(S) # Debugging
    # matprint(L) # Debugging
    # Traversing
    aseq1=[]
    aseq2=[]
    i = [i for i in range(n+1) if alignscore in S[i]][0]
    j = S[i].index(alignscore) 
    while j > 0 and i > 0:
        if S[i][j] == 0 and alignscore !=0:
            break
        # print(L[i][j]) # Debugging
        if L[i][j] == 'd':
            aseq1.append(seq1[i-1])
            aseq2.append(seq2[j-1])
            i-=1
            j-=1
        elif L[i][j] == 'l':
            aseq1.append('-')
            aseq2.append(seq2[j-1])
            j-=1
        elif L[i][j] == 'u':
            aseq1.append(seq1[i-1])
            aseq2.append('-')    
            i-=1
        else:
            i-=1
            j-=1
    # print(aseq1,aseq2) # Debugging
    aseq1.reverse()
    aseq1=''.join(aseq1)
    aseq2.reverse()    
    aseq2=''.join(aseq2)
    return alignscore, aseq1,aseq2

def affine_alignment(seq1,seq2,gapop,gapext,local=False):
    # Needleman-Wunsh & Smith-Whitterman Affine Alignment Algorithms
    # Initialize tables
    n = len(seq1)
    m = len(seq2)  
    G = zeros(n+1,m+1)
    E = zeros(n+1,m+1)
    F = zeros(n+1,m+1)
    V = zeros(n+1,m+1)
    L = zeros(n+1,m+1)
    if local == False:
        # Initialize according to global alignment
        l = -float("inf")
        G[0][0] = -float("inf")
        for j in range(m+1):
            F[0][j] = -float("inf")
        for i in range(n+1):
            E[i][0] = -float("inf")
        for j in range(1,m+1):
            V[0][j] = gapop+j*gapext
        for i in range(1,n+1):
            V[i][0] = gapop+i*gapext
    else:
        # Initialize according to local alignment
        l = 0
        
    for i in range(1,n+1):
        for j in range(1,m+1):
            a = E[i][j-1]+gapext
            b = G[i][j-1]+gapop+gapext
            c = F[i][j-1]+gapop+gapext
            E[i][j]=max(a,b,c)
            
            a = F[i-1][j]+gapext
            b = G[i-1][j]+gapop+gapext
            c = E[i-1][j]+gapop+gapext
            F[i][j]=max(a,b,c)
            
            G[i][j] = V[i-1][j-1]+naive_score(seq1[i-1],seq2[j-1],gapop)
            V[i][j] = max(G[i][j],E[i][j],F[i][j])
            L[i][j] = argmax(G[i][j],E[i][j],F[i][j],l)
    matprint(V) # Debugging
    matprint(L) # Debugging
    if local == False:
        # Score for global alignment
        alignscore = V[-1][-1]
    else:
        # Score for local alignment
         alignscore = max([max(row) for row in V])
    # Traversing
    aseq1=[]
    aseq2=[]
    i = [i for i in range(n+1) if alignscore in V[i]][0]
    j = V[i].index(alignscore)    
    while j > 0 and i > 0:
        # print(L[i][j]) # Debugging
        if V[i][j] == 0 and alignscore !=0:
            break
        if L[i][j] == 'd':
            aseq1.append(seq1[i-1])
            aseq2.append(seq2[j-1])
            i-=1
            j-=1
        elif L[i][j] == 'l':
            aseq1.append('-')
            aseq2.append(seq2[j-1])
            j-=1
        elif L[i][j] == 'u':
            aseq1.append(seq1[i-1])
            aseq2.append('-')    
            i-=1
        else:
            i-=1
            j-=1
    # print(aseq1,aseq2) # Debugging
    aseq1.reverse()
    aseq1=''.join(aseq1)
    aseq2.reverse()    
    aseq2=''.join(aseq2) 
    return alignscore, aseq1,aseq2
"""
# Test & Debug
n1,seq1,n2,seq2 = readfa2("hw3example.fa")
s,aseq1,aseq2 = naive_alignment(seq1,seq2,-3,local=False)
done = writefa2(s,n1,aseq1,n2,aseq2,"hw3example.aln")
"""
# Main Function to Run From Terminal
# allalign --mode alocal --input sequences.fasta --gapopen -5 --gapext -2
try:
    e = 6
    mode =  sys.argv[2] # "--search"
    text_file =  sys.argv[4] #"hw2example.fa"
    gapop = int(sys.argv[6])
    a_n = mode[0]
    if a_n == 'a':
        e = 8
        gapext = int(sys.argv[8])
    else:
        gapext = gapop
except IndexError:    
    print("Insufficient number of arguments! Expected:",e, " Passed:",len(sys.argv)-1)

n1,seq1,n2,seq2 = readfa2(text_file)

if a_n == "a": #Affine?
    try:
        if mode == "aglobal":
            s, aseq1,aseq2 = affine_alignment(seq1,seq2,gapop,gapext,local=False)
            file_name = "global-affineGap.aln"
        elif mode == "alocal":
            s, aseq1,aseq2 = affine_alignment(seq1,seq2,gapop,gapextlocal=True)
            file_name = "local-affineGap.aln"
        else:
            raise NameError
    except NameError:
         print("Unknown Mode:", mode, "\nTry global, local,  aglobal or alocal instead.")
else: #Naive?
    try:
        if mode == "global":
            s,aseq1,aseq2 = naive_alignment(seq1,seq2,gapop,local=False)
            file_name = "global-naiveGap.aln"
        elif mode == "local":
            s,aseq1,aseq2 = naive_alignment(seq1,seq2,gapop,local=True)
            file_name = "local-naiveGap.aln"
        else:
            raise NameError
    except NameError:
         print("Unknown Mode:", mode, "\nTry global, local,  aglobal or alocal instead.")

Done = writefa2(s,n1,aseq1,n2,aseq2,file_name)