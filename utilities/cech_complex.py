#!/usr/bin/env python3
#
# This file is part of 'Aleph - A Library for Exploring Persistent
# Homology'. It permits the calculation of a Čech complex for sets
# of points in the plane.
#
# The input data set is assumed to follow a simple format:
#
#   x_0 x_1 ... x_d
#   y_0 y_1 ... y_d
#
# The previous example contains two points (x,y). Each value needs
# to be separated by whitespace.
#
# Original author: Bastian Rieck

import itertools
import math
import miniball
import sys

if __name__ == "__main__":
  filename  = sys.argv[1]
  radius    = float(sys.argv[2])
  dimension = int(sys.argv[3])
  points    = list()

  with open(filename) as f:
    for line in f:
      coordinates = [ float(x) for x in line.split() ]
      points.append( coordinates )

  n       = len(points)
  indices = list(range(0,n))

  for d in range(0,dimension+1):
    for simplex in itertools.combinations(indices, d):
      if len(simplex) >= 1:
        P = [ tuple( points[i] ) for i in simplex ] # points
        B = miniball.Miniball(P)                    # minimum enclosing ball
        r = math.sqrt( B.squared_radius() )         # radius

        if r <= radius:
          print(simplex)
