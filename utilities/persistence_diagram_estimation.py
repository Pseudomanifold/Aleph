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

def make_prior_ranges(alpha, beta, n, N, k=3):
  mean    = alpha/beta
  var     = alpha/beta**2
  se_mean = mean / math.sqrt(n)
  se_var  = math.sqrt(2*var**2 / (n-1))

  means   = numpy.linspace(mean - se_mean * k, mean + se_mean * k, N)
  vars    = numpy.linspace(var - se_var * k, var + se_var * k, N)

  betas    = [ mean / var for mean,sigma in zip(means,vars) ]
  alphas   = [ mean**2 / var for mean,sigma in zip(means,vars) ]

  return alphas, betas

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

  c_alpha, c_beta = point_estimates( [ c for c,_ in data ] )
  d_alpha, d_beta = point_estimates( [ d for _,d in data ] )

  # Plot the distributions to ensure that the point estimates are
  # actually useful.

  if False:
    v = [ c for c,_ in data ]
    x = numpy.linspace(min(v), max(v), 100)

    plt.hist([c for c,_ in data], normed=True, bins=20, label="Creation [samples]")
    plt.plot(x, pdf(x, d_alpha, c_beta), label="Creation [estimate]")
    plt.legend()
    plt.show()

    v = [ d for _,d in data ]
    x = numpy.linspace(min(v), max(v), 100)

    plt.hist([c for c,_ in data], normed=True, bins=20, label="Destruction [samples]")
    plt.plot(x, pdf(x, d_alpha, d_beta), label="Destruction [estimate]")
    plt.legend()
    plt.show()

  c_alphas, c_betas = make_prior_ranges(c_alpha, c_beta, len(data), 100)

  #v = [ c for c,_ in data ]
  #x = numpy.linspace(min(v), max(v), 100)

  #plt.hist([c for c,_ in data], normed=True, bins=20)

  #for c_alpha, c_beta in zip(c_alphas, c_betas):
  #  plt.plot(x, pdf(x, c_alpha, c_beta))

  #plt.show()
