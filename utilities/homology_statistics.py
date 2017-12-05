#!/usr/bin/env python3
#
# This file is part of 'Aleph - A Library for Exploring Persistent
# Homology'. It contains code for calculating statistics about the
# Betti numbers of a triangulated space.
#
# Original author: Bastian Rieck

import json
import sys

filename = sys.argv[1]

with open(filename) as f:
  data = json.load(f)

homology         = data['homology']
betti_statistics = dict()
euler_statistics = dict()

for index in homology:
  betti = homology[index]['betti'][1]
  euler = homology[index]['euler']

  betti_statistics[betti] = betti_statistics.get(betti, 0) + 1
  euler_statistics[euler] = euler_statistics.get(euler, 0) + 1

print("Betti number statistics:")

for betti in sorted(betti_statistics.keys()):
  print("  %d: %03d" % (betti, betti_statistics[betti]))

print("Euler characteristic statistics:")

for euler in sorted(euler_statistics.keys()):
  print("  %d: %02d" % (euler, euler_statistics[euler]))
