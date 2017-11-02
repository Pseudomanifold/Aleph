#!/usr/bin/env python3
#
# This file is part of the utilities shipped with 'Aleph - A Library for
# Exploring Persistent Homology'. It contains test code to parametrize a
# 'pinched torus'.
#
# Original author: Bastian Rieck

from math import cos, sin, fmod, pi, sqrt

import numpy   as np
import pandas  as pd
import seaborn as sns

import matplotlib.pyplot as plt
import scipy.stats       as ss

n = 4096
m = int(sqrt(n))
R = 10
r = 1

# Gap size in angular coordinates. This is to be seen as the radius for
# which the 'pinch' is relevant.
gap_size = pi / 180.0 * 90

X = list()
Y = list()
Z = list()

for i in range(m):
  for j in range(m):
    phi   = 2*pi * i / (m-1)
    theta = 2*pi * j / (m-1)

    # Angular distance to gap. The tube radius of the torus is modified
    # if this distance is smaller than the gap size.
    dist  = abs(fmod(2*pi+phi,2*pi) - pi)

    if dist <= gap_size:
      # Use a simple linear model to decrease the radius as we are
      # coming towards the gap.
      r_= r / gap_size * dist
    else:
      r_ = r

    x = (R+r_*cos(theta))*cos(phi)
    y = (R+r_*cos(theta))*sin(phi)
    z =  r_ * sin(theta)

    X.append(x)
    Y.append(y)
    Z.append(z)

df = pd.DataFrame(
  {
   'x' : X,
   'y' : Y,
   'z' : Z
  }
)

xmin = min(X)
xmax = max(X)
ymin = min(Y)
ymax = max(Y)
zmin = min(Z)
zmax = max(Z)

kernel = ss.gaussian_kde(np.concatenate([X,Y,Z]),bw_method='silverman')
values = kernel(np.concatenate([X,Y,Z]))

for x,y,z,w in zip(X,Y,Z,values):
  print(x,y,z,w)

sns.pairplot(df, diag_kind="kde")

plt.show()
