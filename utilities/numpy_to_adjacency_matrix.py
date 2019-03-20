#!/usr/bin/env python3
#
# numpy_to_adjacency_matrix.py: converts an `.npy` file to a (textual!)
# adjacency matrix.


import argparse
import sys

import numpy as np


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('INPUT', type=str, help='Input file')
    parser.add_argument('-o', '--output',
        type=str,
        help='Output file (optional)'
    )

    args = parser.parse_args()

    M = np.load(args.INPUT)

    if not args.output or args.output == '-':
        output = sys.stdout
    else:
        output = args.output

    # TODO: make output format configurable?
    np.savetxt(output, M)
