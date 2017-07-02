[![Build Status](https://travis-ci.org/Submanifold/Aleph.svg?branch=master)](https://travis-ci.org/Submanifold/Aleph)

![Aleph logo](Aleph.png "The logo of Aleph in all its glory")

# Aleph &mdash; *A* *L*ibrary for *E*xploring *P*ersistent *H*omology

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

- Clique Community Persistence: A Topological Visual Analysis Approach for Complex Networks
- Hierarchies and Ranks for Persistence Pairs

If you want to contribute, please see the [contribution
guidelines](CONTRIBUTING.md) for more details.

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

* [`FLANN`](https://github.com/mariusmuja/flann) for fast nearest-neighbour queries
* Qt 5, for the GUI

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
you can run the unit test with:

    $ make test

Please submit any issues you may encounter.
