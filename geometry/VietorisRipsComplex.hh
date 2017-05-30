#ifndef ALEPH_GEOMETRY_VIETORIS_RIPS_COMPLEX_HH__
#define ALEPH_GEOMETRY_VIETORIS_RIPS_COMPLEX_HH__

#include "RipsExpander.hh"
#include "RipsSkeleton.hh"

#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include "topology/filtrations/Data.hh"

namespace aleph
{

namespace geometry
{

/**
  Convenience function for building a Vietoris--Rips complex from unstructured
  data. This requires being able to calculate nearest neighbours, as well as a
  maximum connectivity threshold for the complex, and a maximum dimension.

  Note that this function does not give you the opportunity to try out another
  class than the default simplicial complex and simplex class.

  The resulting simplicial complex will use the standard weight function, i.e.
  a simplex has a weight that is equal to the maximum weight of its faces. The
  0-simplices have a weight of 0. 1-simplices use the distance between the two
  vertices of the corresponding edge as a weight.

  Hence, this complex fully represents the scale of the distance function.
*/

template <class NearestNeighbours> auto buildVietorisRipsComplex(
  const NearestNeighbours& nn,
  typename NearestNeighbours::ElementType epsilon,
  unsigned dimension ) -> topology::SimplicialComplex< topology::Simplex<typename NearestNeighbours::ElementType, typename NearestNeighbours::IndexType> >
{
  using ElementType       = typename NearestNeighbours::ElementType;
  using IndexType         = typename NearestNeighbours::IndexType;
  using Simplex           = topology::Simplex<ElementType, IndexType>;
  using SimplicialComplex = topology::SimplicialComplex<Simplex>;

  geometry::RipsSkeleton<NearestNeighbours> ripsSkeleton;

  auto skeleton
    = ripsSkeleton( nn, epsilon );

  geometry::RipsExpander<SimplicialComplex> ripsExpander;

  auto K = ripsExpander( skeleton, dimension );
  K      = ripsExpander.assignMaximumWeight( K );

  K.sort( filtrations::Data<Simplex>() );

  return K;
}

/**
  Convenience function for building a Vietoris--Rips complex from data
  with additional values.

  The additional values are assumed to be specified as a generic range
  of input iterators. The order of values is assumed to be the same as
  the order of points in the wrapper for nearest neighbours.
*/

template <
  class NearestNeighbours,
  class InputIterator
>
auto buildVietorisRipsComplex(
  const NearestNeighbours& nn,
  typename NearestNeighbours::ElementType epsilon,
  unsigned dimension,
  InputIterator begin, InputIterator end ) -> topology::SimplicialComplex< topology::Simplex<typename NearestNeighbours::ElementType, typename NearestNeighbours::IndexType> >
{
  using ElementType       = typename NearestNeighbours::ElementType;
  using IndexType         = typename NearestNeighbours::IndexType;
  using Simplex           = topology::Simplex<ElementType, IndexType>;
  using SimplicialComplex = topology::SimplicialComplex<Simplex>;

  geometry::RipsSkeleton<NearestNeighbours> ripsSkeleton;

  auto skeleton
    = ripsSkeleton( nn, epsilon );

  geometry::RipsExpander<SimplicialComplex> ripsExpander;

  auto K = ripsExpander( skeleton, dimension );
  K      = ripsExpander.assignMaximumData( K, begin, end );

  K.sort( filtrations::Data<Simplex>() );

  return K;
}

} // namespace geometry

} // namespace aleph

#endif
