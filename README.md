![Aleph logo](Aleph.png "The logo of Aleph in all its glory")

# Aleph &mdash; *A* *L*ibrary for *E*ploring *P*ersistent *H*omology

Aleph is a C++ library for exploring and extending usages of [persistent
homology](https://en.wikipedia.org/wiki/Persistent_homology). Its main
goal is to provide users with a versatile, simple-to-use implementation
that quickly permits building prototype applications.

# License

Aleph uses the MIT license. Please see the file [`LICENSE.md`](LICENSE.md)
in the main directory of the repository for more details.

# Requirements

* A recent C++ compiler with support for C++11
* Several `Boost` dependencies for some of the data structures:
  * `Boost.Functional`
  * `Boost.Iterator`
  * `Boost.MultiIndex` 
* [`FLANN`](https://github.com/mariusmuja/flann) for nearest-neighbour queries
