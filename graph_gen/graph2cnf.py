#!/usr/bin/env python

import sys,re

######################################
# OBTAINED FROM https://raw.githubusercontent.com/ibipul/coloring_SAT/master/graph2cnf.py
######################################

# File: To Generate the CNF from the given Graph defination
# Author: Bipul Islam
# Date: May 2013
#usage: ./zchaff graphi.txt>graphi.cnf


for line in open(sys.argv[1]):
    #print line
    
    if line[0]=='p':
        data = line.split()
        colors = int(data[1])
        vertices = int(data[2])
        edges = int(data[3])
        print "p cnf " + str(colors*vertices) + " " + str(vertices + vertices* colors*(colors -1)/2 + colors*edges)
        X = [[0]*(colors+1) for _ in range(vertices+1)]
        counter=1
        i=1
        while( i <vertices+1):
            j=1
            while( j <colors+1):
                          
                X[i][j] = counter
                #print counter
                counter = counter +1        
                    
                j = j + 1
                
            i = i + 1
        
        #print X    
        #print "No of Vertices:: " + vertices + ", No of edges:: " + edges
        i=1;
        while( i<= vertices):
            #print i,
            c=1;
            while( c<= colors):
               # print i,
               # print c,   
               # print "," ,
                print X[i][c],
                if (c==colors):
                    print "0"
                
                    
                c = c + 1
                    
            i = i + 1
        print  
    
        i=1
        while (i <= vertices):
            c=1
            while ( c <= colors-1):
                    d = c + 1
                    while( d <= colors):
                        print "-"+str(X[i][c])+" -"+str(X[i][d])+" 0"
                        d = d+ 1
                
                    c = c + 1
            
            i = i + 1
    
    
    if line[0]=='e':
        data = line.split()
        u = data[1]
        v = data[2]
        
        c=1 
        while ( c <= colors):
            print "-"+str(X[int(u)][c])+ " -" +str(X[int(v)][c])+" 0"
            c = c + 1
