#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys

def setTrans(matrix,From,to,p):
    if(matrix.get(From,None) is None): matrix[From]={}
    matrix[From][to]=p

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(sys.argv[0]+" <in.hmm>\n")
infile=sys.argv[1]
transitions={}
IN=open(infile,"rt")
for line in IN:
    if(line.rstrip()=="transitions:"):
        line=IN.readline()
        fields=line.rstrip().split()
        if(len(fields)!=2): exit("syntax error in hmm file")
        nrow=int(fields[0])
        for i in range(nrow):
            line=IN.readline()
            fields=line.rstrip().split()
            if(len(fields)!=nrow): exit("error in transition matrix")
            for j in range(nrow):
                p=float(fields[j])
                if(p==0): continue
                setTrans(transitions,j,i,p)
IN.close()
keys=list(transitions.keys())
keys.sort(key=lambda x: int(x))
for From in keys:
    column=transitions[From]
    keys2=list(column.keys())
    keys2.sort(key=lambda x: int(x))
    for to in keys2:
        print(str(From)+" -> "+str(to)+" = "+str(column[to]))



