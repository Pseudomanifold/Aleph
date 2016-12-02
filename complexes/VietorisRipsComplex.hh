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
*/

template <class NearestNeighbours> auto buildVietorisRipsComplex(
  const NearestNeighbours& nn,
  typename NearestNeighbours::ElementType epsilon,
  unsigned dimension ) -> SimplicialComplex< Simplex<typename NearestNeighbours::ElementType, typename NearestNeighbours::IndexType> >
{
  using ElementType       = typename NearestNeighbours::ElementType;
  using IndexType         = typename NearestNeighbours::IndexType;
  using Simplex           = Simplex<Element, IndexType>;
  using SimplicialComplex = SimplicialComplex<Simplex>;

  RipsSkeleton<NearestNeighbours> ripsSkeleton;

  auto skeleton
    = ripsSkeleton.build( nn, epsilon );

  RipsExpander<SimplicialComplex> ripsExpander;

  auto K = ripsExpander( skeleton, dimension );
  K      = ripsExpander.assignMaximumWeight( K );

  return K;
}

}

}

#endif
