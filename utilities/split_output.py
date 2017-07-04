#!/usr/bin/env python3
#
# This file is part of 'Aleph - A Library for Exploring Persistent
# Homology'. Its purpose is to split the output of some tools into
# multiple parts, using an appropriate suffix.
#
# Example:
#
# The `network_analysis` tool separates all persistence diagrams by two
# new lines. Say you calculate the persistent homology of a graph:
#
#     $ network_analysis Graph.gml 5 > /tmp/Graph.txt
#
# The file `Graph.txt` will now contain 5 persistence diagrams. To split
# them up, you can call this script as follows:
#
#     $ split_output --prefix=d --digits=1 /tmp/Graph.txt
#
# This will result in the following files:
#
#  - Graph_d0.txt
#  - Graph_d1.txt
#  - Graph_d2.txt
#  - Graph_d3.txt
#  - Graph_d4.txt
#
# You can use any string as a prefix and specify the number of digits to
# use for the prefix generation. By default, two digits will be used, as
# this covers most use-cases. Moreover, you can specify a start value as
# an index using '--start'.
#
# Original author: Bastian Rieck

import argparse
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument('FILE')
parser.add_argument('--prefix', type=str, default="")
parser.add_argument('--digits', type=int, default=2)
parser.add_argument('--start',  type=int, default=0)

arguments = parser.parse_args()
filename  = arguments.FILE

class States:
  Regular, FirstNewLine, SecondNewLine = range(3)

# Keeps track of the current part of the file
index = 0

prefix = arguments.prefix
digits = arguments.digits
start  = arguments.start

def make_output_filename(basename, prefix, digits, index, start):
  return   basename                                  \
         + "_"                                       \
         + prefix                                    \
         + ("%0" + str(digits) + "d") % (index+start)\
         + ".txt"

with open(filename) as f:
  basename = os.path.basename(filename)
  basename = os.path.splitext(basename)[0]

  g = open(make_output_filename(basename, prefix, digits, index, start), "w")

  state = States.Regular
  for line in f:
    if line == "\n":
      if state == States.Regular:
        state = States.FirstNewLine
      elif state == States.FirstNewLine:
        state = States.SecondNewLine

    if state != States.SecondNewLine:
      g.write(line)
    else:
      g.close()
      index = index + 1
      g     = open(make_output_filename(basename, prefix, digits, index, start), "w")
      state = States.Regular

  g.close()
