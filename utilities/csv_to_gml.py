#!/usr/bin/env python3
#
# This file is part of 'Aleph - A Library for Exploring Persistent
# Homology'. It contains code to convert CSV files into graphs, by
# representing them as GML files.
#
# Note that this assumes that the CSV files model a complete graph
# and contain matching rows and columns.
#
# Original author: Bastian Rieck

import csv
import sys

# How many fields to skip in the header of the file. This is used to
# denote the value (name) of the corresponding node.
header_skip_fields = 3

# How many fields to skip in every row. This is used to denote which
# fields contain the name of the corresponding node and the weights,
# which in turn will be used for the edges.
row_skip_fields    = 2

if __name__ == "__main__":

  # Edges are stored in triple form: (u,v,w), where u is the ID of the
  # source edge, v is the ID of the target edge, and w is the weight.
  edges = list()

  with open(sys.argv[1]) as f:
    reader = csv.reader(f)
    header = next(reader)
    names  = header[header_skip_fields:]

    # Map names to indices. This is required because the data file may
    # not necessarily be symmetrical. There may be edges that are only
    # specified in one direction, so they will only occur in the rows.

    last_index    = 0
    name_to_index = dict()

    for index,name in enumerate(names):
      name_to_index[name] = index

    last_index = len(name_to_index)

    # Parse rows. The index is first looked up in the map and updated
    # if necessary.

    for row in reader:
      fields  = row[row_skip_fields:]
      name    = fields[0]
      weights = [ float(x) for x in fields[1:] ]

      if name not in name_to_index:
        name_to_index[name]  = last_index
        last_index          +=1

      u_id   = name_to_index[name]
      u_name = name

      for index,weight in enumerate(weights):
        # The name of the *target* of this edge is guranteed to
        # correspond to item found in the list of names because
        # those names were collected from the columns.
        v_id   = index
        v_name = names[index]

        if u_id != v_id and weight > 0:
          edges.append( (u_id,v_id,weight) )

  ######################################################################
  # Create output
  ######################################################################

  print("graph [")

  for name in sorted(name_to_index.keys()):
    print("  node [\n"
          "    id %d\n"
          "    label \"%s\"\n"
          "  ]" % (name_to_index[name], name))

  for u,v,w in edges:
    print("  edge [\n"
          "    source %d\n"
          "    target %d\n"
          "    weight %.16f\n"
          "  ]" % (u,v,w))

  print("]")
