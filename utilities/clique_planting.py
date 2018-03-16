#!/usr/bin/env python3
#
# This file is part of 'Aleph - A Library for Exploring Persistent
# Homology'. It permits generating a set of graphs for the *clique
# planting problem*.
#
# Original author: Bastian Rieck

import argparse
import random
import sys

import networkx as nx

def plant_clique(G, k):
  """
  Plants a $k$-clique in the graph. To this end, a random subset of the
  vertices is selected and made into a complete subgraph.
  """

  vertices = random.sample( list( range(n) ), k )

  for index, u in enumerate(vertices):
    for v in vertices[index+1:]:
      G.add_edge(u,v)

  return G

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Generate clique planting problem instances.')

  parser.add_argument('--n', type=int, default=1000, help='Number of vertices in the graph')
  parser.add_argument('--k', type=int, default=  25, help='Size of planted clique'         )

  arguments = parser.parse_args()

  n = arguments.n
  k = arguments.k
  p = 0.5
  G = nx.erdos_renyi_graph(n, p)

  # With probability = 0.5, we decide upon whether a clique should be
  # planted in the graph or not.
  if random.choice( [True, False] ):
    G = plant_clique(G, k)

  nx.write_gml(G, sys.stdout.buffer)
