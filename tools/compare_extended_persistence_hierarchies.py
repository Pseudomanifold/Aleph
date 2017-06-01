#!/usr/bin/env python3

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
