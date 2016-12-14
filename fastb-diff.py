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
from Fastb import Fastb

if(len(sys.argv)!=4):
    exit(sys.argv[0]+" <from.fastb> <minus.fastb> <out.fastb>")
(fromFile,minusFile,outFile)=sys.argv[1:]

fromFastb=Fastb(fromFile)
minusFastb=Fastb(minusFile)
n=fromFastb.numTracks()
for i in range(n):
    fromTrack=fromFastb.getIthTrack(i)
    id=fromTrack.getID()
    minusTrack=minusFastb.getTrackByName(id)
    fromData=fromTrack.getData()
    minusData=minusTrack.getData()
    L=len(fromData)
    if(len(minusData)!=L): exit("track lengths differ")
    for j in range (L): fromData[j]-=minusData[j]
fromFastb.save(outFile)


