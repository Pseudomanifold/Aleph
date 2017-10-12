#!/usr/bin/env python3
#
# simple_network_analysis.py: Performs simple methods for network
# analysis in order to establish a reasonable "baseline" for more
# sophisticated methods.
#
# This file is part of 'Aleph - A Library for Exploring Persistent
# Homology'. It is meant for stand-alone usage.
#
# Original author: Bastian Rieck

import numpy    as np
import networkx as nx
import pandas   as pd

import os
import sys

if __name__ == "__main__":
  n         = len(sys.argv)-1
  filenames = sys.argv[1:n-1]
  labels    = np.loadtxt(sys.argv[n])
  graphs    = list()

  for filename in filenames:
    if os.path.splitext(filename)[1] == ".gml":
      G = nx.read_gml(filename, label='id')
      graphs.append(G)

  ######################################################################
  # Associate features with the graph
  ######################################################################
  #
  # - Average clustering coefficient
  # - Average degree of nodes
  # - Average shortest path length
  # - Density
  # - Diameter

  df = pd.DataFrame(columns=['average_clustering_coefficient',\
                             'average_degree',\
                             'average_shortest_path_length',\
                             'density',\
                             'diameter'])

  for G in graphs:
    average_clustering_coefficient = nx.average_clustering(G)
    average_degree                 = np.mean( [degree for _,degree in nx.degree(G) ] )

    average_shortest_path_length = 0
    diameter                     = 0

    for g in nx.connected_component_subgraphs(G):
      average_shortest_path_length = max(average_shortest_path_length, nx.average_shortest_path_length(g))
      diameter                     = max(diameter, nx.diameter(g))

    density = nx.density(G)
    data    = {
      'average_clustering_coefficient': average_clustering_coefficient,
      'average_degree'                : average_degree,
      'average_shortest_path_length'  : average_shortest_path_length,
      'density'                       : density,
      'diameter'                      : diameter
    }

    df = df.append(data, ignore_index=True)
