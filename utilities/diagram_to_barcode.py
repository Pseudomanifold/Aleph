#!/usr/bin/env python3
#
# This file is part of 'Aleph - A Library for Exploring Persistent
# Homology'. It contains code for converting a persistence diagram
# into a simple persistence barcode.
#
# The format follows common plotting programs such as `gnuplot` or
# `pgfplots` (for LaTeX) in order to enable rapid usage in papers,
# technical reports, or other documents.

import sys

data = list()

with open(sys.argv[1]) as f:
  for line in f:
    line  = line.strip()

    if not line:
      continue

    (x,y) = [ float(t) for t in line.split() ]
    data.append( (x,y) )

data = sorted(data)

for i,(x,y) in enumerate(data):
  x0,y0 = x,i
  x1,y1 = y,i

  print("{}\t{}\n{}\t{}\n".format(x0,y0,x1,y1))
