[![Build Status](https://travis-ci.org/Submanifold/Aleph.svg?branch=master)](https://travis-ci.org/Submanifold/Aleph) [![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/972/badge)](https://bestpractices.coreinfrastructure.org/projects/972)

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

## Getting started with Aleph

Aleph is a header-only library. Nonetheless, it ships with numerous
tests, examples, and tools that require compilation. Coincidentally,
these files also demonstrate how to use the library in practice.

See the following links for more details:

- [Building & installing Aleph](building.md): how to build the library,
  its tests and associated tools, and install it afterwards
- [Publications](publications.md): see how Aleph was used in some
  publications and how to reproduce some results
- [Examples](examples.md): well-documented & small example programs that
  demonstrate various functions
- [Tutorials](tutorials.md): how to solve certain common tasks in
  topological data analysis with Aleph
