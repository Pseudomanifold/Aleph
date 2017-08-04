#!/usr/bin/env python3
#
# Extracts documentation from a set of directories and stores them as
# reStructuredText files. The reason for this script is that I prefer
# using Sphinx for the documentation but the C++ interface is not yet
# capable of documenting an API. Rather than *separating* source from
# documentation, I want to go a middle way here.
#
# Note that this script does *not* perform any kind of syntax checks,
# but uses the extracted comments "as-is" and concatenates them.

import os
import re
import sys

from textwrap import dedent

# TODO: make configurable?
directories = [ "../include", "../src" ]
extensions  = [ ".cc", ".hh" ]
output_root = "../docs/sphinx/source/API"

re_docstring = re.compile(r'/\*\*(?:[^*]|\*(?!/))*\*/', re.MULTILINE)

"""
Extracts documentation from the given filename and puts a file with the
same name into the root output directory.
"""
def extract_documentation(root, filename, output_root):
  output = ""
  data   = ""

  with open( os.path.join(root, filename ) ) as f:
    data = f.read()

  for s in re_docstring.findall(data):
    # Prune the result prior to writing it to the file: remove the
    # starting portion of the comment and the end portion. This is
    # easy because we just have to drop a few characters.
    s = s[4:]
    s = s[:-2]
    s = dedent(s)

    output += s

  if output:
    output_filename = os.path.splitext(filename)[0] + ".rst"
    with open(os.path.join(output_root, output_filename), "w") as f:
      f.write(output)

for directory in directories:
  for root, dirs, files in os.walk(directory):
    for filename in files:
      if os.path.splitext(filename)[1] in extensions:
        extract_documentation(root, filename, output_root)
