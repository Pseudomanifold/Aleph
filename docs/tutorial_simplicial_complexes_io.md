---
title: Reading and writing simplicial complexes
---

# Reading and writing simplicial complexes

Usually, simplicial complexes are not built manually by
users&nbsp;(Aleph offers this option, however, in case you need it;
please see [the corresponding tutorial](tutorial_simplicial_complex_manually.md)
for more information). Instead, they may be read from a variety of different
file formats, such as certain *graph file formats* and *meshes*. If you
want to build a simplicial complex from point cloud data, for example
via Vietoris&ndash;Rips expansion, please take a look at
[this tutorial](tutorial_vietoris_rips_complex_point_cloud.md) instead.

# Supported file formats

Aleph uses individual reader classes for every supported file format.
Please take a look at `include/aleph/topology/io` for more information.
Currently, the following formats are supported:

- Edge lists, the simplest and most common file format for exchanging
  graph data. Aleph supports reading optional *weights* and *node
  labels*.

- [`HDF5`](https://support.hdfgroup.org/HDF5), the *Hierarchical Data
  Format*. Aleph supports reading simple *data spaces* from HDF5 files.

- [`GML`](https://gephi.org/users/supported-graph-formats/gml-format),
  the *graph modeling language* format. Notice that Aleph supports an
  extremely liberal superset of this language and is much more forgiving
  in parsing certain dialects than other software such as [NetworkX](https://networkx.github.io).

- [`Pajek`](http://vlado.fmf.uni-lj.si/pub/networks/pajek), also known
  as the *NET graph format*. Aleph supports reading&nbsp;(weighted)
  graphs in this format.

- [`PLY`](http://paulbourke.net/dataformats/ply), the *Stanford Polygon
  File Format*. Aleph supports reading an ASCII variant of this file
  format, which is commonly used to describe 3D meshes.

- [`VTK`](https://www.vtk.org/doc/nightly/html/classvtkStructuredGrid.html),
  the file format introduced by the *Visualization Toolkit*. At present,
  Aleph supports reading *structured grids*.

In case a reader does not support your needs, please open an issue in
the project's [issue
tracker](https://github.com/Submanifold/Aleph/issues)&mdash;I am always
interested in extending Aleph's capabilities.

# The easy way

The easy way to load a simplicial complex from a supported file format
involves using the class `aleph::topology::io::SimplicialComplexReader`.
This class is a wrapper for *all* supported formats and provides an
easy interface:

```cpp
#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/io/SimplicialComplexReader.hh>

// Let's specify a simplicial complex first. If you want to *guess* the
// file format, the process is a little bit more involved.
using DataType          = double;
using VertexType        = unsigned;
using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

// [...]


// The reader class has to know where data is to be stored. Errors will
// be raised depending on some readers (some formats do not offer this,
// because their limited specification does not permit error handling).
SimplicialComplex K;

// The correct reader is subsequently chosen depending on the extension
// of the input file.
auto filename = "test.gml";

aleph::topology::io::SimplicialComplexReader reader;
reader( filename, K );
```

Additional information *may* be available&mdash;depending on the data
format, of course&mdash;and can be queried:

```cpp
// The indexing of the labels follows the vertex order of the simplicial
// complex. This works even if the ordering does not begin with zero.
auto labels = reader.labels();

for( auto&& label : label )
  std::cout << label << "\n";
```

Please refer to the documentation of this class for more information.
Notice that the interface of the reader class is by necessity somewhat
limited because it cannot provide all possible options of all readers.

# Reading a specific file format

If you are only interested in loading a specific file format, it is
easiest to use the corresponding reader class from
`include/aleph/topology/io`. For example, suppose you want to load
graphs in the common *edge lists* format, i.e. a specification in which
edges are stored as two vertex indices, followed by an optional weight.
The `EdgeListReader` class is suitable for this purpose and offers some
customization options:

```cpp
#include <aleph/topology/io/EdgeLists.hh>

// With the same `using` declarations as above...

aleph::topology::io::EdgeListReader reader;
reader.setReadWeights( true ); // we expect edge weights

reader.setSeparator( ":" );    // we expect vertex indices to be separated
                               // by colons instead of whitespace

reader.setTrimLines( true );   // we want to throw away and ignore all other
                               // information
```

Other readers may offer even more attributes and behavioural changes.
Please refer to the documentation for more details.

# More information

The `SimplicialComplexReader` class is used in various tests:

* [`test_io_gml.cc`](https://github.com/Submanifold/Aleph/blob/master/tests/test_io_gml.cc)
* [`test_io_pajek.cc`](https://github.com/Submanifold/Aleph/blob/master/tests/test_io_pajek.cc)
* [`test_io_vtk.cc`](https://github.com/Submanifold/Aleph/blob/master/tests/test_io_vtk.cc)

A short demonstration is also available in the [`relevant_edges.cc`](https://github.com/Submanifold/Aleph/blob/master/src/tools/relevant_edges.cc) tool.
