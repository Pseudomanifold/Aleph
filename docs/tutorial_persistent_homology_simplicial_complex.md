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

auto persistenceDiagrams = calculatePersistenceDiagram<Standard>( K );

for( auto&& diagram : persistenceDiagrams )
  std::cout << diagram << "\n";
```

This should result in the same persistence diagram. The test suite of
Aleph, in particular [`test_persistent_homology_complete.cc`](https://github.com/Submanifold/Aleph/blob/master/tests/test_persistent_homology_complete.cc) defines numerous tests to ensure that neither the representation nor the reduction algorithm influence the results&nbsp;(with the possible exception of runtime and memory usage, of course).
