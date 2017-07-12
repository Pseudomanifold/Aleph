---
title: Examples
---

Here's a list of examples that are present and documented in Aleph. Each example comes with a link
that leads directly to the corresponding source file, as well as a brief description. A good start
for an absolute novice would be the Vietoris&mdash;Rips calculation example.

| Name | Description |
|------|-------------|
| [create_persistence_diagrams](https://github.com/Submanifold/Aleph/blob/master/examples/create_persistence_diagrams.cc) | Creating random persistence diagrams in numerous ways |
| [create_random_graph](https://github.com/Submanifold/Aleph/blob/master/examples/create_random_graph.cc) | Creating a random graph and storing it as a GML file |
| [network_analysis](https://github.com/Submanifold/Aleph/blob/master/examples/network_analysis.cc) | Analysing the persistent homology of a weighted network (a graph) |
| [ply](https://github.com/Submanifold/Aleph/blob/master/examples/ply.cc) | Calculating persistent homology of meshes in ASCII PLY format with different filtrations; includes statistical information about persistence diagrams |
| [vietoris_rips](https://github.com/Submanifold/Aleph/blob/master/examples/vietoris_rips.cc) | Calculating a Vietoris&mdash;Rips complex of point clouds |
| [vietoris_rips_eccentricity](https://github.com/Submanifold/Aleph/blob/master/examples/vietoris_rips.cc) | Calculating a Vietoris&mdash;Rips complex of point clouds using an eccentricity data descriptor |
| [vtk](https://github.com/Submanifold/Aleph/blob/master/examples/vtk.cc) | Calculating persistent homology of structured grids in VTK legacy format |

# Building the examples

To build the examples, please follow the [instructions on how to build
Aleph](building.md). Examples will automatically be built and placed in
the `examples` subfolder of the build directory.
