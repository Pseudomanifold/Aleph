#!/usr/bin/env python3
#
# This file is part of 'Aleph - A Library for Exploring Persistent
# Homology'. It converts a matrix defined in an implicit format as
# used by `gnuplot`, for example, into a surface with coordinates,
# so that it can be used with additional software like TikZ.
#
# All output is written to STDOUT.
#
# Original author: Bastian Rieck

import re
import sys

rows = 0
cols = 0

with open(sys.argv[1]) as f:
  for line in f:
    line = line.strip()
    cols = max(cols, len( line.split() ))
    rows = rows+1

print("Identified (%d,%d)-matrix" % (rows, cols), file=sys.stderr)

with open(sys.argv[1]) as f:
  row = 0
  for line in f:
    line   = line.strip()
    values = line.split()
    
    for col in range(cols):
      print("%d %d %s" % (col,row,values[col]))

    print("")
    row = row+1

print("Finished conversion", file=sys.stderr)
