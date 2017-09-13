#!/usr/bin/env python3
#
# This file is part of the utilities shipped with 'Aleph - A Library for
# Exploring Persistent Homology'. Its purpose is the conversion of image
# data to matrix data.
#
# Usage: image_to_matrix FILE
#
# The matrix will be written to `stdout` in matrix format. Note that the
# function does *not* perform any conversions. Hence, it works best with
# simple greyscale images.
#
# Original author: Bastian Rieck

import numpy
import sys

from PIL import  Image

image = Image.open(sys.argv[1])
array = numpy.array(image)

numpy.savetxt(sys.stdout.buffer, array, fmt="%s")
