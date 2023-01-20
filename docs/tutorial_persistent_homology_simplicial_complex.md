---
title: Calculating persistent homology of a simplicial complex
---

# Calculating persistent homology of a simplicial complex

Previously, we learned [how to build simplicial complexes
manually](tutorial_simplicial_complex_manually.md). We will now see how
to calculate the persistent homology of such a complex.

## Direct calculation

In the following, we assume that we are given a simplicial complex `K`
whose persistent homology we want to calculate. If simplices already
have weights that are consistent with respect to the filtration
ordering, we can immediately obtain some results:

```cpp
#include <aleph/persistentHomology/Calculation.hh>

auto persistenceDiagrams = calculatePersistenceDiagram( K );

for( auto&& diagram : persistenceDiagrams )
  std::cout << diagram << "\n";
```

In this example, the current ordering of the simplicial complex is being
used as the filtration. The calculation of persistent homology uses the
default algorithms, as specified in `aleph/config/Defaults.hh`. The file
may look like this, for example:

```cpp
#include <aleph/persistentHomology/algorithms/Twist.hh>

#include <aleph/topology/representations/Vector.hh>

namespace aleph
{

namespace defaults
{

using Index              = unsigned;
using Representation     = topology::representations::Vector<Index>;
using ReductionAlgorithm = persistentHomology::algorithms::Twist;

} // namespace defaults

} // namespace aleph
```

We can see that, by default, the *twist* algorithm&nbsp;(see
[*Persistent homology computation with
a twist*](https://eurocg11.inf.ethz.ch/abstracts/22.pdf) for more
details) is used to reduce the boundary matrix of the simplicial
complex.

In its technical innards, Aleph always uses a boundary matrix
representation for calculating persistent homology. This matrix may be
represented in various forms. Here, an `std::vector` is used to
store one column of the matrix. Aleph borrows the concept of multiple
representations from [*PHAT&mdash;Persistent Homology Algorithms Toolbox*](https://people.mpi-inf.mpg.de/~mkerber/bkrw-pphat.pdf). Normally, you do not have to change this.

However, it may make sense to change the reduction algorithm. To use the
standard persistent homology calculation algorithm, use the following
call:

```cpp
using namespace aleph::persistentHomology::algorithms;

auto persistenceDiagrams = calculatePersistenceDiagrams<Standard>( K );

for( auto&& diagram : persistenceDiagrams )
  std::cout << diagram << "\n";
```

This should result in the same persistence diagram. The test suite of
Aleph, in particular [`test_persistent_homology_complete.cc`](https://github.com/Pseudomanifold/Aleph/blob/master/tests/test_persistent_homology_complete.cc) defines numerous tests to ensure that neither the representation nor the reduction algorithm influence the results&nbsp;(with the possible exception of runtime and memory usage, of course).

## Filtrations & expansion

The case described above rarely occurs in practice. Often, we first need
to assign weights to a simplicial complex, expand it, and bring it into
filtration order. Let us return to the previous data set, a simple
triangle, for which we only specify vertices and edges. We then expand
the resulting simplicial complex to a Vietoris&ndash;Rips complex. This
will automatically create the triangle for us:

```cpp
#include <aleph/topology/Simplex.hh>

using namespace aleph::topology;

using Data              = double;
using Vertex            = unsigned;
using Simplex           = Simplex<Data, Vertex>;
using SimplicialComplex = SimplicialComplex<Simplex>;

std::vector<Simplex> simplices
  = { {1}, {2}, {4}, {1,2}, {1,4}, {2,4} };

SimplicialComplex K( simplices.begin(), simplices.end() );
RipsExpander<SimplicialComplex> ripsExpander;

auto L = ripsExpander( K, 2 );

std::cout << L.size()              << "\n"; // 7; the last simplex is the triangle
std::cout << L.contains( {4,2,1} ) << "\n"; // true
```

The simplicial complex `L` now contains the triangle `{4,2,1}` because
the simplicial complex `K` already contains all of its edges. At this
point, we do *not* have any weights attached to the complex, though, so
let us add some:

```cpp
// We could use any sequence type here, of course. The weights have to
// follow the ordering of vertices. Hence, the assignment will be:
//
// vertex -> data
// 1      -> 1
// 2      -> 2
// 4      -> 3
std::vector<DataType> data = {
  1, 2, 3
};

L = ripsExpander.assignMaximumData( L, data.begin(), data.end() );
```

By calling `assignMaximumData()`, we instruct the expander class to
assign a simplex `s` the maximum data of all of its faces. In total, we
will have the following mapping:

    {1}     -> 1
    {2}     -> 2
    {4}     -> 3
    {1,2}   -> 2
    {1,4}   -> 3
    {2,4}   -> 3
    {1,2,4} -> 3

Next, we need to sort the simplicial complex according to these weights:

```cpp
L.sort( aleph::topology::filtrations::Data<Simplex>() );
```

This is a standard sorting predicate that sorts simplices by their data
values and, in case of ties, lexicographically. Hence, lower-dimensional
simplices are guaranteed to preceded higher-dimensional ones. Having
established the filtration order, we may proceed to calculate persistent
homology as usual:

```cpp
auto persistenceDiagrams = calculatePersistenceDiagrams( L );

for( auto&& diagram : persistenceDiagrams )
  std::cout << diagram << "\n";
```

Note that `assignMaximumData()` is a convenience function for assigning,
well, the maximum data value to simplices. It corresponds to a sublevel
set filtration of the simplicial complex. Other filtrations are of
course possible. This is managed by a flexible functor-based system. To
use it, we have to provide an initial value that is to be used when
determining the data value of a simplex, and a functor. To obtain
a superlevel set filtration, you can use the following code, for
example:

```cpp
L = ripsExpander.assignData( L, data.begin(), data.end(),
                             std::numeric_limits<Data>::max(),
                             [] ( const Data& a, const Data& b )
                             {
                               return std::min(a,b);
                             } );
```

Thanks to the power of lambda expressions in C++11, you can use even
more complex data assignments for the simplicial complex. Just make sure
that the resulting filtration is consistent&mdash;faces *must* precede
co-faces here!

## More information

Please take a look at the test suite of Aleph for more information.
Persistent homology calculation is well-tested for different scenarios.
The following tests may prove useful:

- [`test_connected_components.cc`](https://github.com/Pseudomanifold/Aleph/blob/master/tests/test_connected_components.cc)
- [`test_persistent_homology_complex.cc`](https://github.com/Pseudomanifold/Aleph/blob/master/tests/test_persistent_homology_complete.cc)
- [`test_rips_expansion.cc`](https://github.com/Pseudomanifold/Aleph/blob/master/tests/test_rips_expansion.cc)

In the [next tutorial](tutorial_vietoris_rips_complex_point_cloud.md), we will put everything together and cover the
complete Vietoris&ndash;Rips complex expansion process of a point cloud.
