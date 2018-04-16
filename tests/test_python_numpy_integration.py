#!/usr/bin/env python3

import aleph as al

try:
    import numpy as np
except ImportError:
    return

import unittest

class TestSimplexMethods(unittest.TestCase):
    def test_construction(self):
        pass

M = al.SimplicialComplex([[0], [1], [2], [1, 0], [2, 0]])

diagram = al.calculatePersistenceDiagrams(M)[0]
numpy_diagram = np.array(diagram)

for point, np_point in zip(diagram, numpy_diagram):
    assert point.x == np_point[0]
    assert point.y == np_point[1]
