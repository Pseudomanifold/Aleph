#!/usr/bin/env python3

from aleph import *

import unittest

class TestSimplexMethods(unittest.TestCase):
  def test_construction(self):
    pass

a   = Simplex( [0] )
b   = Simplex( [1] )
c   = Simplex( [2] )
ab  = Simplex( [0,1] )
ac  = Simplex( [0,2] )
bc  = Simplex( [1,2] )
abc = Simplex( [0,1,2] )

K = SimplicialComplex( [a,b,c,ab,ac,bc,abc] )
L = SimplicialComplex( [ [0], [1], [0,1] ] )
M = SimplicialComplex( [ ([0], 0), ([1], 2), ([2], 0) ] )

print(K)
print(L)
print(M)

for vertex in abc:
  print(vertex)

for simplex in abc.boundary:
  print(simplex)

assert 0 in abc    , "Vertex not found"
assert 3 not in abc, "Incorrect vertex"

print( len(abc) )
print( abc == abc )
print( a.dimension )

print( len(K) )
print( len(L) )
print( K.dimension )
print( L.dimension )
print( "First simplex:", K[0] )
print( "Last simplex:", K[ len(K) -1 ] )

print( "Does the complex contain the simplex?", abc in K )
print( "Does the complex not contain the simplex?", Simplex( [4711] ) not in K )

if abc:
  print("Non-empty simplex")

print( abc[2] )

diagrams = calculatePersistenceDiagrams(K)
for diagram in diagrams:
  print(diagram)

print(sorted(M, key=lambda x: x.data ))

def sorting_function(s,t):
  return s.data < t.data

M.sort( sorting_function )

print(M)

norms.pNorm(diagrams[-1])
norms.totalPersistence(diagrams[-1])
norms.infinityNorm(diagrams[-1])
