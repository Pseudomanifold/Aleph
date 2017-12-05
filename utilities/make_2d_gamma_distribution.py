#!/usr/bin/env python3
#
# This file is part of 'Aleph - A Library for Exploring Persistent
# Homology'. It creates a set of points in the plane, whose values
# follow two different Gamma distributions.
#
# Original author: Bastian Rieck

import numpy as np

if __name__ == "__main__":
  n      = 20
  alphas = [7.5,9.0]
  betas  = [1.0,2.0]

  for i in range(n):
    for j in range(n):
      x = np.random.gamma(alphas[0], 1.0 / betas[0])
      y = np.random.gamma(alphas[1], 1.0 / betas[1])

      print(x,y)
