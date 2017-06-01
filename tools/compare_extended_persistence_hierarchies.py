#!/usr/bin/env python3

import numpy
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
          if destroyer != float('inf'):
            destroyer = int(destroyer * factor)
          else:
            destroyer = 1e16

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

filenames   = sys.argv[1:]
hierarchies = list()

for filename in filenames:
  hierarchies.append( load_hierarchy(filename, scale=True, factor=1000) )

n = len(hierarchies)
M = numpy.zeros( (n,n) )

for i in range(n):
  for j in range(i+1,n):
    d = zss.simple_distance(hierarchies[i], hierarchies[j],
                            PersistencePair.get_children,
                            PersistencePair.get_label,
                            PersistencePair.distance)

    M[i,j] = d
    M[j,i] = d

numpy.savetxt("/tmp/M.txt", M)
