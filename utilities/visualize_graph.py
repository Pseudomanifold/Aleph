#!/usr/bin/env python3
#
# This file is part of 'Aleph - A Library for Exploring Persistent
# Homology'. Its purpose is to draw a set of simple graphs using a
# circular layout. The graphs are supposed to be in GML format and
# the script is particularly geared towards visualizing discussion
# graphs in which a single node is highly important.
#
# The script will accept a number of input graphs, and store their
# visualizations in '/tmp'.
#
# Original author: Bastian Rieck

import networkx          as nx
import matplotlib.pyplot as plt

import os
import sys

filenames = sys.argv[1:]
for filename in filenames:
  G  = nx.read_gml( filename, label='id' )
  Gc = max(nx.connected_component_subgraphs(G), key=len)
  A  = nx.nx_agraph.to_agraph(Gc)

  print("Processing %s..." % filename)

  # Get the root node for the subsequent layout
  D  = dict(nx.degree(Gc))
  r  = max(D.keys(), key=lambda x: D[x])

  A.graph_attr['root'] = r
  A.node_attr['shape'] = 'point'

  A.layout(prog='twopi')

  name = os.path.basename(filename)
  name = os.path.splitext(name)[0]

  A.draw("/tmp/"+name+".svg")
