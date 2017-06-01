#!/usr/bin/env python3

import re
import sys
import zss

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

  # Returns the 'label' of the node, which merely consists of a tuple
  # with the creator and destroyer information.
  @staticmethod
  def get_label(node):
    return (node.creator, node.destroyer)

  # Calculates the distance between two nodes.
  @staticmethod
  def distance(u, v):

    if not u:
      return abs(v[0]-v[1])
    elif not v:
      return abs(u[0]-u[1])

    c1,c2 = u[0],v[0]
    d1,d2 = u[1],v[1]

    return max(abs(c1-c2),abs(d1-d2))

"""
Reads an extended persistence hierarchy from a file and converts it to
a sequence of nodes, as specified above. During the conversion process
the values may optionally be scaled. If 'scale' is set and 'factor' is
valid, then the 'raw' floating point values read from the file will be
multiplied by 'factor' and afterwards converted into an int. This step
is required because the zss module currently only supports ints during
the distance calculations.
"""
def load_hierarchy(filename, scale=False, factor=None):
  reNode = r'(\d+):\s+(\d+)\s+(\S+)'
  reEdge = r'(\d+)\s+--\s+(\d+)'

  id_2_pair = dict()
  edges     = list()

  with open(filename) as f:
    for line in f:
      match = re.match(reNode, line)
      if match:
        id        = int(match.group(1))
        creator   = float(match.group(2))
        destroyer = float(match.group(3))

        if scale and factor:
          creator   = int(creator * factor)
          destroyer = int(destroyer * factor)

        id_2_pair[id] = PersistencePair(creator, destroyer)
      else:
        match = re.match(reEdge, line)
        if match:
          edges.append( (int(match.group(1)),int(match.group(2))) )

A = PersistencePair(6,100).add_child(
    PersistencePair(4,3)
  ).add_child(
    PersistencePair(5,2)
  )

B = PersistencePair(6,100).add_child(
    PersistencePair(4,3).add_child(
      PersistencePair(5,2)
    )
    )

print(zss.simple_distance(A,B, PersistencePair.get_children, PersistencePair.get_label, PersistencePair.distance))
