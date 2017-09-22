---
title: Building Aleph
---

# Building Aleph

Aleph uses the [CMake](https://cmake.org) build system. In the best
case, all dependencies will automatically be identified. Please use
your package manager to install the following dependencies:

* A recent C++ compiler with support for C++11 and, preferably, C++14.
  Aleph is being continuously built by `g++`, the C++ compiler of the
  GNU project, and `clang++`, the C++ compiler of the LLVM project.
  Other compiler may or may not work. Please open an
  [issue](https://github.com/Submanifold/Aleph/issues) if you require
  another compiler.
* `CMake`, preferably a recent version >= 3.2
* Several `Boost` dependencies for some of the data structures:
  * `Boost.Functional`
  * `Boost.Iterator`
  * `Boost.MultiIndex`

Aleph is meant to be used as a header-only library on top of which you
can develop your own projects based on persistent homology. However,
Aleph ships with numerous unit tests, some example programs, and tools
required for my current research. For building them, please clone the
repository to some local directory on your computer. Running the
following commands within this directory should be sufficient in most
cases:

    $ mkdir build
    $ cd build
    $ cmake ../
    $ make

It is advisable to test that Aleph works correctly on your system, so
you can run the unit tests with:

    $ make test

A typical output of this command should look like this:

    Running tests...
    Test project /home/bastian/Projects/Aleph/build
          Start  1: bootstrap
     1/22 Test  #1: bootstrap ........................   Passed    0.00 sec
          Start  2: boundary_matrix_reduction
     2/22 Test  #2: boundary_matrix_reduction ........   Passed    0.00 sec
          Start  3: clique_enumeration
     3/22 Test  #3: clique_enumeration ...............   Passed    0.00 sec
          Start  4: clique_graph
     4/22 Test  #4: clique_graph .....................   Passed    0.00 sec
          Start  5: connected_components
     5/22 Test  #5: connected_components .............   Passed    0.05 sec
          Start  6: data_descriptors
     6/22 Test  #6: data_descriptors .................   Passed    0.23 sec
          Start  7: filesystem
     7/22 Test  #7: filesystem .......................   Passed    0.00 sec
          Start  8: graph_generation
     8/22 Test  #8: graph_generation .................   Passed    0.00 sec
          Start  9: io_functions
     9/22 Test  #9: io_functions .....................   Passed    0.00 sec
          Start 10: io_gml
    10/22 Test #10: io_gml ...........................   Passed    0.02 sec
          Start 11: io_pajek
    11/22 Test #11: io_pajek .........................   Passed    0.01 sec
          Start 12: io_vtk
    12/22 Test #12: io_vtk ...........................   Passed    0.63 sec
          Start 13: mesh
    13/22 Test #13: mesh .............................   Passed    0.00 sec
          Start 14: munkres
    14/22 Test #14: munkres ..........................   Passed    0.00 sec
          Start 15: nearest_neighbours
    15/22 Test #15: nearest_neighbours ...............   Passed    0.10 sec
          Start 16: persistence_diagrams
    16/22 Test #16: persistence_diagrams .............   Passed    5.50 sec
          Start 17: persistent_homology_complete
    17/22 Test #17: persistent_homology_complete .....   Passed    0.18 sec
          Start 18: point_clouds
    18/22 Test #18: point_clouds .....................   Passed    0.09 sec
          Start 19: rips_expansion
    19/22 Test #19: rips_expansion ...................   Passed    0.00 sec
          Start 20: rips_skeleton
    20/22 Test #20: rips_skeleton ....................   Passed    0.04 sec
          Start 21: step_function
    21/22 Test #21: step_function ....................   Passed    0.00 sec
          Start 22: union_find
    22/22 Test #22: union_find .......................   Passed    0.00 sec

    100% tests passed, 0 tests failed out of 22

    Total Test time (real) =   6.87 sec

Please submit any issues you may encounter. If you want to help, please
take a look at the [contribution guidelines](https://github.com/Submanifold/Aleph/blob/master/CONTRIBUTING.md).

# Additional options

Some of the components of Aleph may be disabled if you want to increase
performance of the build. Currently, these components include the
examples and the tools. Use `ccmake .` in the build directory to get an
overview of the available options. A typical output of this command
would be:

    BUILD_EXAMPLES        ON
    BUILD_PYTHON_BINDINGS ON
    BUILD_TOOLS           ON

Toggle the required options and re-configure the project in order to
change what is built.

# Installing Aleph

Having built Aleph, you can install the library as follows:

    $ make install

It is possible to change the installation prefix of Aleph using either
`cmake` or `ccmake`, the console front-end. To set another prefix upon
building Aleph, you may use the following command:

    $ cmake -DCMAKE_INSTALL_PREFIX=/path/to/installation/directory ../
    $ make install

Package files (for certain Linux distributions) are forthcoming. Please
refer to [this issue](https://github.com/Submanifold/Aleph/issues/27)
for more details.
