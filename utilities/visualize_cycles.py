#!/usr/bin/env python3
#
# This file is part of 'Aleph - A Library for Exploring Persistent
# Homology'. It visualizes information about simple cycles in some
# set of graphs, providing a sort of mean histogram visualization,
# as the time or simulation parameter of the graph varies.
#
# Original author: Bastian Rieck

import matplotlib.cm     as cm
import matplotlib.pyplot as plt
import numpy             as np

import csv
import sys

class Entry:
  def __init__(self, t, num_cycles, cycle_lengths):
    self.t             = t
    self.num_cycles    = num_cycles
    self.cycle_lengths = cycle_lengths

  def __repr__(self):
    output = "t:             {}\n"\
             "num_cycles:    {}\n"\
             "cycle_lengths: {}\n".format(self.t, self.num_cycles, " ".join( str(length) for length in self.cycle_lengths))

    return output

if __name__ == "__main__":
  filename   = sys.argv[1]
  t_min      =   sys.maxsize  # first time step
  t_max      =  -sys.maxsize  # last time step
  min_length =   sys.maxsize  # minimum cycle length
  max_length =  -sys.maxsize  # maximum cycle length

  with open(filename) as f:
    reader = csv.reader(f)

    # Stores individual entries from the CSV file along with their
    # properties such as the number of cycles.
    data = dict()
    for row in reader:
      t             = int(row[0])
      num_cycles    = int(row[1])
      cycle_lengths = [int(length) for length in row[2].split()]
      t_min         = min(t, t_min)
      t_max         = max(t, t_max)
      min_length    = min(min_length, min(cycle_lengths))
      max_length    = max(max_length, max(cycle_lengths))

      data[t] = Entry(t, num_cycles, cycle_lengths)

  num_time_steps = len(data)
  num_rows       = max_length - min_length + 1
  num_cols       = num_time_steps
  data_array     = np.zeros((num_rows, num_cols), dtype=float)

  for index,t in enumerate(sorted(data.keys())):
    cycle_lengths = data[t].cycle_lengths
    for length in cycle_lengths:
      row_index = length - min_length
      col_index = index

      data_array[row_index, col_index] += 1

  data_array = ( data_array - np.min(data_array, axis=0) ) / np.max(data_array, axis=0)

  plt.yticks(np.arange(0, max_length-min_length+1,1), np.arange(min_length, max_length+1,1))
  plt.imshow(data_array, origin="lower", aspect=10, cmap=cm.magma)
  plt.show()
