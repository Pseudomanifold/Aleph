#ifndef ALEPH_COMPLEXES_RIPS_SKELETON_HH__
#define ALEPH_COMPLEXES_RIPS_SKELETON_HH__

#include "Simplex.hh"
#include "SimplicialComplex.hh"

#include <vector>

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
    auto numVertices = nn.size();

    std::vector<Simplex> simplices;
    simplices.reserve( numVertices );

    for( decltype(numVertices) i = 0; i < numVertices; i++ )
      simplices.push_back( Simplex( i ) );

    return SimplicialComplex( simplices.begin(), simplices.end() );
  };
};

}

}

#endif
