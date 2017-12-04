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
  filenames  = sys.argv[1:]
  t_min      =   sys.maxsize  # first time step
  t_max      =  -sys.maxsize  # last time step
  min_length =   sys.maxsize  # minimum cycle length
  max_length =  -sys.maxsize  # maximum cycle length

  # Global data, aggregated over *all* files. No normalization will be
  # performed when performing the aggregation per file.
  data = dict()

  for filename in filenames:
    print("Processing file '{}'...".format(filename), file=sys.stderr)
    with open(filename) as f:
      reader = csv.reader(f)

      # Stores individual entries from the CSV file along with their
      # properties such as the number of cycles.
      local_data = dict()
      for row in reader:
        t             = int(row[0])
        num_cycles    = int(row[1])
        cycle_lengths = [int(length) for length in row[2].split()]
        t_min         = min(t, t_min)
        t_max         = max(t, t_max)
        if cycle_lengths:
          min_length    = min(min_length, min(cycle_lengths))
          max_length    = max(max_length, max(cycle_lengths))

        local_data[t] = Entry(t, num_cycles, cycle_lengths)

        # Aggregate local data and global data
        if t not in data.keys():
          data[t] = local_data[t]
        else:
          num_cycles    = data[t].num_cycles + local_data[t].num_cycles
          cycle_lengths = data[t].cycle_lengths + local_data[t].cycle_lengths
          data[t]       = Entry(t, num_cycles, cycle_lengths)

  num_time_steps = len(data)
  num_rows       = max_length - min_length + 1
  num_cols       = num_time_steps
  data_array     = np.zeros((num_rows, num_cols), dtype=float)

  print("Maximum cycle length: {}\n"
        "Minimum cycle length: {}\n".format(max_length, min_length), file=sys.stderr)

  for index,t in enumerate(sorted(data.keys())):
    cycle_lengths = data[t].cycle_lengths
    for length in cycle_lengths:
      row_index = length - min_length
      col_index = index

      data_array[row_index, col_index] += 1

  data_array = ( data_array - np.min(data_array, axis=0) ) / np.max(data_array, axis=0)
  data_array = np.nan_to_num(data_array) # Replace NaNs with zero. They can only
                                         # be generated in cases where the `max`
                                         # function returns zero, so we may just
                                         # as well replace them.

  # Output to STDOUT; this could be made configurable with respect to
  # the formatting specifier.
  np.savetxt(sys.stdout.buffer, data_array, fmt="%.5f")

  fig, axes = plt.subplots()

  plt.yticks(np.arange(0, max_length-min_length+1,1), np.arange(min_length, max_length+1,1))
  plt.imshow(data_array, origin="lower", aspect=5, cmap=cm.magma)
  plt.colorbar()
  plt.show()
