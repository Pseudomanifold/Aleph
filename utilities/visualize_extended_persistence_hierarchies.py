#!/usr/bin/env python3
#
# This file is part of 'Aleph - A Library for Exploring Persistent
# Homology'. It contains code for visualizing extended persistence
# hierarchies (also called interlevel set persistence hierarchies)
# as TikZ pictures.
#
# The file processes _all_ command-line arguments and expects them
# to be hierarchies. The generated code will be written to STDOUT.
#
# For more information, please refer to:
#
#    Hierarchies and Ranks for Persistence Pairs
#    Bastian Rieck, Heike Leitte, and Filip Sadlo
#    Proceedings of TopoInVis 2017, Japan

import re
import sys

"""
Class describing a persistence pair in the extended persistence
hierarchy. This class is basically a container for all critical
points of a function.
"""
class PersistencePair:
  def __init__(self, creator, destroyer):
    self.creator   = creator
    self.destroyer = destroyer
    self.children  = list()

  # Adds a new child to the current node; this function return 'self' in
  # order to permit chaining (not that we use it in this example).
  def add_child(self, child):
    self.children.append(child)
    return self

  @staticmethod
  def get_children(node):
    return node.children

"""
Reads an extended persistence hierarchy from a file and converts it to
a sequence of nodes, as specified above.
"""
def load_hierarchy(filename, scale=False, factor=None):
  reNode = r'(\d+):\s+([\d\.\d]+)\s+(\S+)'
  reEdge = r'(\d+)\s+--\s+(\d+)'

  id_2_pair = dict()
  edges     = list()

  with open(filename) as f:
    for line in f:
      match = re.match(reNode, line)
      if match:
        id        = int(match.group(1))
        creator   = int(match.group(2))
        destroyer = float(match.group(3))

        if destroyer == float('inf'):
          destroyer = int(1e4) # TODO: hard-coded...

        id_2_pair[id] = PersistencePair(creator, destroyer)
      else:
        match = re.match(reEdge, line)
        if match:
          edges.append( (int(match.group(1)),int(match.group(2))) )

  for id_u,id_v in edges:
    u = id_2_pair[id_u]
    v = id_2_pair[id_v]

    u.add_child(v)

  # Find root(s): a root is a node that only appears as the source of an
  # edge but not as the destination
  sources = set([u for u,_ in edges])
  targets = set([v for _,v in edges])
  roots   = sources.difference(targets)

  assert len(roots) == 1, "Hierarchy must be connected"
  return id_2_pair[roots.pop()]

"""
Formats a persistence pair for TikZ output.
"""
def format_pair(pair):
  return "node { $(%d,%d)$ }" % (pair.creator, pair.destroyer)

"""
Traverses and prints the hierarchy.
"""
def traverse_hierarchy(node, level):
  prefix = ' ' * level * 3
  output = ''

  for index,child in enumerate(node.children):
    output += prefix
    output += "child {\n"
    output += prefix + '  '
    output += format_pair(child)
    output += "\n"
    output += traverse_hierarchy(child, level+1)
    output += prefix + '  ' + "}"
    if index+1 == len(node.children) and level == 1:
        output += ";\n"
    else:
      output += "\n"

  return output

filenames   = sys.argv[1:]
hierarchies = list()

for filename in filenames:
  hierarchies.append( load_hierarchy(filename) )

for hierarchy in hierarchies:
  print("\\begin{tikzpicture}")
  print("  \\%s" % format_pair(hierarchy) )
  print("    %s" % traverse_hierarchy(hierarchy, 1))
  print("\\end{tikzpicture}")
