#ifndef ALEPH_GEOMETRY_VIETORIS_RIPS_COMPLEX_HH__
#define ALEPH_GEOMETRY_VIETORIS_RIPS_COMPLEX_HH__

#include "Simplex.hh"
#include "SimplicialComplex.hh"

#include "RipsExpander.hh"
#include "RipsSkeleton.hh"

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
  unsigned dimension ) -> SimplicialComplex< Simplex<typename NearestNeighbours::ElementType, typename NearestNeighbours::IndexType> >
{
  using ElementType       = typename NearestNeighbours::ElementType;
  using IndexType         = typename NearestNeighbours::IndexType;
  using Simplex           = Simplex<ElementType, IndexType>;
  using SimplicialComplex = SimplicialComplex<Simplex>;

  complexes::RipsSkeleton<NearestNeighbours> ripsSkeleton;

  auto skeleton
    = ripsSkeleton.build( nn, epsilon );

  complexes::RipsExpander<SimplicialComplex> ripsExpander;

  auto K = ripsExpander( skeleton, dimension );
  K      = ripsExpander.assignMaximumWeight( K );

  return K;
}

}

}

#endif
