#!/usr/bin/env python3
#
# This file is part of 'Aleph - A Library for Exploring Persistent
# Homology'. It counts elementary cycles in a graph and prints the
# count to STDOUT. The script tries to be smart and uses *numbers*
# detected in the input filenames in order to structure its output
# accordingly.
#
# Original author: Bastian Rieck

import networkx as nx

import csv
import os
import re
import sys

if __name__ == "__main__":
  filenames = sys.argv[1:]

  for filename in filenames:
    G       = nx.read_gml(filename, label="id")
    basis   = nx.cycle_basis(G)
    matches = re.match(r'[^0-9]+(\d+).*', os.path.basename(filename))

    if matches:
      t             = matches.group(1)
      num_cycles    = len(basis)
      cycle_lengths = [len(cycle) for cycle in basis]
      row           = [num_cycles]
      row.extend(cycle_lengths)

      # Output to STDOUT; this is more flexible and permits us to use
      # information in a variety of scripts.
      writer = csv.writer(sys.stdout)
      writer.writerow(row)
