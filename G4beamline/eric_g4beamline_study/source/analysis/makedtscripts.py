#!/usr/bin/python
# This file will generate a series of script files to be used by g4beamline
# usage  python makescripts.py [baseFile] [minDelta] [maxDelta]
import os
import sys

#default values
minDelta = 2.0
maxDelta = 10.00
baseFile = sys.argv[1]                        # base file name
if (len(sys.argv)>2):
  minDelta = float(sys.argv[2])
if (len(sys.argv)>3):
  maxDelta = float(sys.argv[3])+.00001        # can get into trouble with roundoff

subsString = "xxxx"                           # string to substitute

deltas = []
delta = minDelta
while (delta<=maxDelta):
  deltaString = "%.1f"%delta
  outputFile = baseFile.replace(subsString,deltaString)
  os.system("sed 's/%s/%s/g' %s > %s"%(subsString,deltaString,baseFile,outputFile))
  if (delta<3.0):
    delta+=0.1
  else:
    delta+=1.
    

