#ifndef ALEPH_COMPLEXES_RIPS_SKELETON_HH__
#define ALEPH_COMPLEXES_RIPS_SKELETON_HH__

#include "Simplex.hh"
#include "SimplicialComplex.hh"

namespace aleph
{

namespace complexes
{

template <class NearestNeighbours> class RipsSkeleton
{
public:
  using ElementType       = typename NearestNeighbours::ElementType;
  using IndexType         = typename NearestNeighbours::IndexType;

  using Simplex           = Simplex<ElementType, IndexType>;
  using SimplicialComplex = SimplicialComplex<Simplex>;

  SimplicialComplex build( const NearestNeighbours& nn, ElementType epsilon ) const
  {
  };
};

}

}

#endif
