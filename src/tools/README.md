This directory contains some tools that are built using the `Aleph`
framework. While you should be able to compile all of them ideally,
the `Makefile` exports individual targets. Hence, every tool can be
built separately.

# Requirements

To build most of the tools, you should install at least the following
packages:

- A modern C++ compiler that supports (at least) the C++11 standard
- The `Boost` C++ libraries
- `CMake`

# Building a specific tool

To build a specific tool, clone this repository and execute the
following commands in the cloned directory:

```console
$ mkdir build
$ cd build
$ cmake ../
$ make TOOL_NAME
```

Here `TOOL_NAME` is one of the internal tool names, as specified in the
headers of the subsequent sections.

# `connectivity_matrix_analysis`

This tool can be used to analyse the persistent homology of connectivity
matrices, stored in text format. To convert such a matrix from a general
`numpy` format, you should use the `numpy_to_adjacency_matrix.py` script
from the `utilities` directory:

```console
python3 numpy_to_adjacency_matrix.py -o matrix.txt matrix.npy
```

This will create a file called `matrix.txt` that can be read by the tool
via the following commands from the `build` directory:

```console
./tools/connectivity_matrix_analysis -i 2 matrix.txt
```

This will create a `JSON` file with topological descriptor data. Here is
an example of such a file (with fake numbers):

```json
{
	"diagrams": [{
		"betti": 0,
		"dimension": 0,
		"persistent_entropy": "6.593758",
		"total_persistence_1": "108.269235",
		"total_persistence_2": "129.286312",
		"name": "matrix.txt",
		"size": 2,
		"diagram": [
			["0", "1.07421"],
			["0", "1.48229"],
		]
	}]
}
```

The entries `persistent_entropy` and `total_persistence_*` are scalar
summaries of the corresponding persistence diagram than could be used
for a first assessment of structural differences.
