#!/usr/bin/env python3

import math
import numpy
import sys

import matplotlib.pyplot as plt
import scipy.stats       as stats

def point_estimates(samples):
  mean  = numpy.mean(samples)
  var   = numpy.var(samples)
  alpha = mean * mean / var
  beta  = mean / var

  return alpha,beta

def pdf(x, alpha, beta):
  return stats.gamma.pdf(x, a=alpha, scale=1.0/beta)

if __name__ == "__main__":
  filename = sys.argv[1]
  data     = list()

  with open(filename) as f:
    for line in f:
      if not line.strip() or line.startswith('#'):
        continue

      (c,d) = [ float(x) for x in line.strip().split() ]

      if not math.isinf(c) and not math.isinf(d):
        data.append( (c,d) )

  alpha_c, beta_c = point_estimates( [ c for c,_ in data ] )
  alpha_d, beta_d = point_estimates( [ d for _,d in data ] )

  # Plot the distributions to ensure that the point estimates are
  # actually useful.

  v = [ c for c,_ in data ]
  x = numpy.linspace(min(v), max(v), 100)

  plt.hist([c for c,_ in data], normed=True, bins=20, label="Creation [samples]")
  plt.plot(x, pdf(x, alpha_c, beta_c), label="Creation [estimate]")
  plt.legend()
  plt.show()

  v = [ d for _,d in data ]
  x = numpy.linspace(min(v), max(v), 100)

  plt.hist([c for c,_ in data], normed=True, bins=20, label="Destruction [samples]")
  plt.plot(x, pdf(x, alpha_d, beta_d), label="Destruction [estimate]")
  plt.legend()
  plt.show()
