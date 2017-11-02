#!/usr/bin/env python3

from math import cos
from math import pi
from math import sin
from math import sqrt

import pandas  as pd
import seaborn as sns

import matplotlib.pyplot as plt

n = 1024
m = int(sqrt(n))
R = 10
r = 1

X = list()
Y = list()
Z = list()

for i in range(m):
  for j in range(m):
    phi   = 2*pi * i / (m-1)
    theta = 2*pi * j / (m-1)

    x = (R+r*cos(theta))*cos(phi)
    y = (R+r*cos(theta))*sin(phi)
    z =  r * sin(theta)

    print(x,y,z)

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

sns.pairplot(df)

plt.show()
