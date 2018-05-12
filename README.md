[![Build Status](https://travis-ci.org/Submanifold/Aleph.svg?branch=master)](https://travis-ci.org/Submanifold/Aleph) [![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/972/badge)](https://bestpractices.coreinfrastructure.org/projects/972)

![Aleph logo](Aleph.png "The logo of Aleph in all its glory")

# Aleph &mdash; A Library for Exploring Persistent Homology

Aleph is a C++ library for exploring and extending usages of [persistent
homology](https://en.wikipedia.org/wiki/Persistent_homology). Its main
goal is to provide users with a versatile, simple-to-use implementation
that quickly permits building prototype applications.

Aleph is inspired by [`DIPHA`](https://github.com/DIPHA/dipha) and
[`PHAT`](https://bitbucket.org/phat-code/phat). In particular, Aleph
borrowed the idea of keeping the *representation* of a boundary matrix
separate from the implementation.

For more information, please read the [original paper describing
`PHAT`](https://people.mpi-inf.mpg.de/~mkerber/bkrw-pphat.pdf).

Since its inception in late 2016, Aleph has been used to support the
following papers:

- [Clique Community Persistence: A Topological Visual Analysis Approach for Complex Networks](https://submanifold.github.io/Aleph/Rieck17d.html)
- [Persistence Concepts for 2D Skeleton Evolution Analysis](https://submanifold.github.io/Aleph/Rieck17b.html)
- [Hierarchies and Ranks for Persistence Pairs](https://submanifold.github.io/Aleph/Rieck17a.html)
- [*Shall I compare thee to a network?*&mdash;Visualizing the Topological
  Structure of Shakespeare&rsquo;s plays](https://submanifold.github.io/Aleph/Rieck16b.html)

Please refer to [the list of
publications](https://submanifold.github.io/Aleph/publications) in the
documentation of Aleph for more details. The documentation covers how to
reproduce a subset of the results reported in the papers above.

If you want to contribute, please see the [contribution
guidelines](CONTRIBUTING.md) for more details.

# Features

Aleph contains numerous algorithms and helper classes that simplify
working with persistent homology. Here is a brief selection of the most
important ones:

* Easy-to-use and expressive simplex and simplicial complex class
* Support for different input formats to read simplicial complexes from
  a variety of input files
    - 1D functions
    - Edge lists
    - GraphML
    - GML
    - HDF5
    - Lexicographic triangulations
    - Matrices
    - NET (Pajek graphs)
    - PLY
    - VTK
* Standard algorithm and *twisted* reduction algorithm for calculating
  persistent homology
* Support for *dualized* variants of both algorithms
* Support for different boundary matrix representations
* Persistence diagram class
* Distance and kernel measures
    - Bottleneck distance
    - Multi-scale smoothing kernel
    - Wasserstein distance
* Algorithms for computing [*intersection homology*](http://www.math.ias.edu/~goresky/pdf/IH.pdf) and [*persistent intersection homology*](https://doi.org/10.1007/s10208-010-9081-1)
* Basic support for *Čech complexes*
* Support for *Dowker complexes*

# Documentation

[Documentation](https://submanifold.github.io/Aleph) of the main
features, including some tutorials, is available on GitHub. If you want
to delve into the code, the `examples` subdirectory is a good starting
point.

# License

Aleph uses the MIT license. Please see the file [`LICENSE.md`](LICENSE.md)
in the main directory of the repository for more details.

# Requirements

* A recent C++ compiler with support for C++11
* `CMake`, preferably a recent version >= 3.2
* Several `Boost` dependencies for some of the data structures:
  * `Boost.Functional`
  * `Boost.Iterator`
  * `Boost.MultiIndex`

# Optional dependencies

* [`Eigen3`](http://eigen.tuxfamily.org) for some auxiliary mathematical
  functions
* [`FLANN`](https://github.com/mariusmuja/flann) for fast nearest-neighbour queries
* [`HDF5`](https://www.hdfgroup.org) for parsing `HDF5` input files
* [`pybind11`](https://github.com/pybind/pybind11) for building the
  Python bindings
* [`RapidJSON`](http://rapidjson.org) for parsing JSON input files
* [`TinyXML2`](https://github.com/leethomason/tinyxml2) for parsing
  GraphML input files

# Building Aleph

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

Please submit any issues you may encounter.

For more information, including how to run tests, please refer to the
[detailed build
instructions](https://submanifold.github.io/Aleph/building) in the
documentation.

# Installing Aleph

If you want to install Aleph from source, simply issue

    $ make install

from the compilation directory. In general, this will require `root`
privileges, unless you change the `CMAKE_INSTALL_PREFIX` variable to
a *local directory*.

It is easier to install Aleph as a package. Currently, only packages
for Arch Linux are available. Use your favourite AUR helper tool for
installing Aleph:

    $ pacaur -S aleph-git # pacaur (deprecated)
    $ trizen -S aleph-git # trizen
    $ yaourt -S aleph-git # yaourt

If you want to volunteer and submit a package for your favourite Linux
distribution, please take a look at [issue #27](https://github.com/Submanifold/Aleph/issues/27)
and add your comments.

# Installing the Python bindings

If your build instructions are configured to build the Python bindings,
i.e. `BUILD_PYTHON_BINDINGS` follow these instructions to install them:

    $ cd build/bindings/python/aleph
    $ pthon3 setup.py install

Note that this uses the old `setuptools` approach for installing the
package. A simpler installation based on `pip` is forthcoming.

# Contact

For general discussion, questions, and comments, please contact the
principal developer and maintainer Bastian Rieck (bastian.rieck@iwr.uni-heidelberg.de).
