#!/usr/bin/env python3
#
# This is an examnple file shipped by 'Aleph - A Library for Exploring
# Persistent Homology'.
#
# This example demonstrates how to load a network, i.e. a graph, from
# a variety of input files. The graph will subsequently be *expanded*
# to a simplicial complex, whose persistent homology is calculated.
#
# This example demonstrates the Python interface of Aleph. Note that the
# functionality presented in this file should be roughly equivalent to:
#
#   examples/network_analysis.cc
#
# Some deviations between the C++ interface and the Python interface are
# allowed and to be expected.
#
# Original author: Bastian Rieck

from aleph import *

import argparse
import sys

if __name__ == "__main__":

  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--invert-weights",    action="store_true", help="if specified, inverts weights stored in the network")
  parser.add_argument("-n", "--normalize-weights", action="store_true", help="if specified, normalizes weights stored in the network to [0,1]")
  parser.add_argument("FILE", type=str, help="input file; must contain a *valid* network that can be read by Aleph")

  arguments = parser.parse_args()
  filename  = arguments.FILE

  sys.stderr.write("Loading '%s'..." % filename)
  sys.stderr.flush()

  K = load(filename)

  sys.stderr.write("finished\n")

  ######################################################################
  # Determine weight range
  ######################################################################

  def get_min_max_weights(K):
    weights = [ s.data() for s in K ]
    return min(weights), max(weights)

  get_min_max_weights(K)

  # TODO:
  #  - Rips expansion
  #  - Weight assignment

  diagrams = calculatePersistenceDiagrams( K )

  for diagram in diagrams:
    diagram.removeDiagonal()

    for point in diagram:
      print(point)

    print(diagram)
