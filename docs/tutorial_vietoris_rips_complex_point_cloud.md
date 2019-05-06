---
title: Calculating a Vietoris&ndash;Rips complex of a point cloud
---

# Calculating a Vietoris&ndash;Rips complex of a point cloud

Building a [Vietoris&ndash;Rips
complex](https://en.wikipedia.org/wiki/Vietoris%E2%80%93Rips_complex) is
a central operation for topological data analysis. In previous
tutorials, we saw how to [build simplicial complexes
manually](tutorial_simplicial_complex_manually.md), [load them from
files](tutorial_simplicial_complexes_io.md), and [calculate their
persistent
homology](tutorial_persistent_homology_simplicial_complex.md).

This tutorial demonstrates how to obtain a simplicial complex from an
unstructured point cloud.

## Point clouds

In the terminology of Aleph, a *point cloud* is any unstructured set of
data that can be represented in tabular form. Aleph follows the
assumption that data points are contained in the *rows* of said table,
while the *columns* contain data attributes. The default representation
of these point clouds is `aleph::containers::PointCloud`. This class is
templated over the type of data to store and expects that each attribute
has the same (numerical) type.

At present, Aleph only supports loading point clouds in ASCII format.
Individual values need to be separated by at least one separator
character. Valid separator characters are:

- any form of whitespace (` `)
- a comma (`,`)
- a semicolon (`;`)
- a colon (`:`)

Hence, Aleph is capable of reading well-formatted numerical CSV files,
for example. Other example point clouds, such as the famous [Iris data
  set](https://en.wikipedia.org/wiki/Iris_flower_data_set) in several
  variants, are provided in the [`tests` directory](https://github.com/Pseudomanifold/Aleph/tree/master/tests/input).

To load a point cloud, use `aleph::containers::load()`:

```cpp
#include <aleph/containers/PointCloud.hh>

// You need to know/specify the attribute type before-hand. This cannot
// be changed while the program is running!
using DataType = double;

auto filename   = "test.csv";
auto pointCloud = aleph::containers::load<DataType>( filename );
```

The point cloud class offers some standard query functions:

```cpp
std::cout << pointCloud.dimension() << "\n"  // dimension of each point
          << pointCloud.size()      << "\n"  // number of points
          << pointCloud.empty()     << "\n"; // does the point cloud contain data?
```

Individual points can be set or stored in various ways:

```cpp
// Let's assume that the point cloud is 3-dimensional. If this is not
// the case, some of the methods below will throw.
pointCloud.set(i, {1,2,3} );

std::vector<DataType> p = {4,5,6};
std::list<DataType>   q = {7,8,9};

pointCloud.set(i+1, p.begin(), p.end() );
pointCloud.set(i+2, q.begin(), q.end() );

// The default get method merely requires an output iterator, making it
// highly-generic.
q.clear();
pointCloud.get(i+1, std::back_inserter(q));

// Element access is also possible. This mimicks the way one accesses
// values in an `std::vector` class.
p = pointCloud[i];
```

## Neighbourhood calculations 

So far, we have only *loaded* the point cloud. The Vietoris&ndash;Rips
calculation necessitates an expansion based on neighbourhoods, though.
To this end, Aleph, offers a set of wrapper classes. The idea is to make
it easy to test out different methods of calculating neighbourhoods,
depending on the dimensionality of the input data, for example. In
high-dimensional spaces, approximations to the nearest neighbours are
often extremely useful, while in low-dimensional spaces, brute-force
calculations are often sufficient.

Let us briefly take a look at two wrappers for neighbourhood
calculations. First, the brute-force wrapper. It has no external
dependencies but is extremely slow for larger data sets and higher
dimensions.

```cpp
#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/distances/Euclidean.hh>

using DataType   = double;
using PointCloud = aleph::containers::PointCloud<DataType>;
using Distance   = aleph::geometry::distances::Euclidean;

aleph::geometry::BruteForce<PointCloud, Distance> bruteForceWrapper;
```

The preceding declaration shows that the wrapper requires
a *container*&nbsp;(in this case, the point cloud) and
a *distance*&nbsp;(in this case, the Euclidean distance). The choice of
distance measure is very important with respect to the structures that
can be detected. Aleph provides some distances measures that permit
calculating distances between two ranges of iterators. Please refer to
the [source code](https://github.com/Pseudomanifold/Aleph/tree/master/include/aleph/geometry/distances) for a list of valid distance measures.

Another useful wrapper is the wrapper based on
[FLANN](https://github.com/mariusmuja/flann), the *Fast Library for
Approximate Nearest Neighbors*. If this library is installed on your
system, Aleph detects it upon installation and makes it automatically
available. For compatibility reasons, it is recommended that you
*always* check for the availability of this library, though:

```cpp
#include <aleph/config/FLANN.hh>

#include <aleph/geometry/FLANN.hh>
#include <aleph/geometry/distances/Euclidean.hh>

using DataType   = double;
using PointCloud = aleph::containers::PointCloud<DataType>;
using Distance   = aleph::geometry::distances::Euclidean;

#ifdef ALEPH_WITH_FLANN 
aleph::geometry::FLANN<PointCloud, Distance> flannWrapper;
#endif
```

If FLANN is available, `ALEPH_WITH_FLANN` will be defined in the
configuration file.

## Vietoris&ndash;Rips expansion

Each wrapper offers the same interface, permitting a *radius search* in
order to obtain all points within a given radius of each other, as well
as a *neighbour search* in order to to obtain the, say, five nearest
neighbours of every point. For the purposes of this tutorial, it is not
required that we use this interface directly&mdash;luckily, Aleph
provides a simple interface for calculating Vietoris&ndash;Rips
complexes of a given radius:

```cpp
#include <aleph/geometry/VietorisRipsComplex.hh>

// some more definitions and `using` declarations as above

aleph::geometry::BruteForce<PointCloud, Distance> bruteForceWrapper;

auto K = aleph::geometry::buildVietorisRipsComplex( bruteForceWrapper,
                                                    epsilon,
                                                    dimension );
```

This convenience function uses the `epsilon` parameter to determine
a maximum search radius with respect to the selected distance measure.
Moreover, it requires a maximum dimension up to which the simplicial
complex is to be expanded. However, the function cannot guarantee that
the expansion contains simplices of the maximum dimension&mdash;this
obviously depends on the data set. A good value for a point cloud is
given by `pointCloud.dimension() + 1`.

The previous code automatically uses the distance measure in order to
define the weight for one-simplices and higher-dimensional simplices,
while the vertices are assigned a weight of zero. Hence, the resulting
simplicial complex is a model of the distance function on the data.

## Using function values for functional persistent homology

Often, one wants to use a set of function values, defined on the points
of the point cloud, to assign to the Vietoris&ndash;Rips complex. For
example, in my paper [*Exploring and Comparing Clusterings of
Multivariate Data Sets Using Persistent
Homology*](http://bastian.rieck.ru/research/EuroVis2016.pdf), I am using
the results of&nbsp;(essentially) a density estimator to assign function
values to the point cloud. Those function values are in turn then used
as weights for the simplicial complex. In his paper [*Topological
pattern recognition for point cloud
data*](https://doi.org/10.1017/S0962492914000051), Gunnar Carlsson
refers to this process as *functional persistence*.

Aleph provides a convenience function for this purpose as well. You need
a range of function values whose cardinality matches the number of
vertices in the simplicial complex. Values in this range will then be
assigned to the vertices of the simplicial complex upon expanding it.
The order of this assignment follows the lexicographical ordering of
vertices in the complex. This means that in virtually all cases, the
function just works exactly the way you expect it to. 

Let us try this for a simple example:

```cpp
// again, let's assume that the `using` declarations from above are
// still present

std::vector<DataType> functionValues = {1,2,3,...};

aleph::geometry::BruteForce<PointCloud, Distance> bruteForceWrapper;

auto K = aleph::geometry::buildVietorisRipsComplex( bruteForceWrapper,
                                                    epsilon,
                                                    dimension,
                                                    functionValues.begin(), functionValues.end() );
```

The simplicial complex will automatically be ordered by increasing
function values. Moreover, each simplex is assigned the maximum function
value of its faces.

## Other weight assignment schemes

You can also assign weights manually. To this end, just iterate over the
simplicial complex and use the `replace()` function to replace
individual simplices.

## Changing the filtration

As mentioned above, the simplicial complex is automatically ordered by
increasing simplex weights. Aleph offers other filtrations in order to
satisfy different topological data analysis scenarios.

An *upper-star filtration*, for example, is obtained as follows:

```cpp
#include <aleph/topology/filtrations/UpperStar.hh>

// The weights of the simplicial complex need to be recalculated
// if we are to use another filtration.
//
// There is a nice convenience function for this purpose. To use
// it, we specify false as the first argument---this indicates a
// re-calculation based on the minimum function value of an edge
// or higher-dimensional simplex.
//
// The second optional parameter of this function can be used to
// skip 1-dimensional simplices (edges).
K.recalculateWeights( false );

// Function values to use for the upper-star filtration below
std::vector<DataType> eccentricity = {1,2,3,...};

// The upper-star filtration is implemented as a sorting
// predicate. It uses the eccentricity values calculated
// above and prepares a lookup table for each simplex.
aleph::topology::filtrations::UpperStar<Simplex> upperStarFiltration( eccentricity.begin(), eccentricity.end() );

// Note that for reasons of simplicity, the sorting predicate is
// always copied during the sorting. To prevent copies, we hence
// use a reference wrapper from the STL.
K.sort( std::ref( upperStarFiltration ) );
```

See the [eccentricity-based Vietoris&ndash;Rips complex calculation
example](https://github.com/Pseudomanifold/Aleph/blob/master/examples/vietoris_rips_eccentricity.cc)
for more information.

## Calculating persistent homology

Given the Vietoris&ndash;Rips complex, the calculation of persistent
homology, depending on the currently-selected filtration, is very
simple:

```cpp
#include <aleph/persistentHomology/Calculation.hh>

auto persistenceDiagrams = aleph::calculatePersistenceDiagrams( K );

for( auto&& diagram : persistenceDiagrams )
{
  // Removes all features of zero persistence. They only clutter up
  // the diagonal.
  diagram.removeDiagonal();

  std::cout << "Persistence diagram:\n"
            << diagram << "\n\n";
}
```

The function `calculatePersistenceDiagrams()` is a convenience function
that uses the default reduction algorithm to calculate persistent
homology. Everything is configurable via templates, though, so you could
for example force the function to use the standard reduction algorithm
by doing the following:

```cpp
#include <aleph/persistentHomology/algorithms/Standard.hh>

using ReductionAlgorithm = aleph::persistentHomology::algorithms::Standard;

auto persistenceDiagrams = aleph::calculatePersistenceDiagrams<Standar>( K );
```

In general, the default choices should be fine for most application
scenarios. Note that the choice of algorithm only affects the
performance but *not* the results! Aleph tests each reduction algorithm
meticulously in order to ensure that the persistence diagrams remain the
same!

## More information

For more information, the following&nbsp;(fully-documented) examples may
prove useful:

- [Vietoris&ndash;Rips expansion with distances](https://github.com/Pseudomanifold/Aleph/blob/master/examples/vietoris_rips.cc)
- [Vietoris&ndash;Rips expansion with eccentricity data descriptor](https://github.com/Pseudomanifold/Aleph/blob/master/examples/vietoris_rips_eccentricity.cc)
