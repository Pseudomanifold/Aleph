#ifndef ALEPH_GEOMETRY_RIPS_SKELETON_HH__
#define ALEPH_GEOMETRY_RIPS_SKELETON_HH__

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <vector>

namespace aleph
{

namespace geometry
{

template <class NearestNeighbours> class RipsSkeleton
{
public:
  using ElementType       = typename NearestNeighbours::ElementType;
  using IndexType         = typename NearestNeighbours::IndexType;

  using Simplex           = topology::Simplex<ElementType, IndexType>;
  using SimplicialComplex = topology::SimplicialComplex<Simplex>;

  SimplicialComplex operator()( const NearestNeighbours& nn, ElementType epsilon ) const
  {
    auto numVertices = nn.size();

    std::vector<Simplex> simplices;
    simplices.reserve( numVertices );

    for( decltype(numVertices) i = 0; i < numVertices; i++ )
      simplices.push_back( Simplex( i ) );

    std::vector< std::vector<IndexType> > indices;
    std::vector< std::vector<ElementType> > distances;

    nn.radiusSearch( epsilon, indices, distances );

    for( std::size_t i = 0; i < indices.size(); i++ )
    {
      IndexType u = IndexType(i);

      for( std::size_t j = 0; j < indices[i].size(); j++ )
      {
        IndexType v   = IndexType( indices[i][j] );
        ElementType d = distances[i][j];

        // Ensures that edges are only added once. Otherwise, edge u--v will
        // also give rise to v--u.
        if( u < v )
          simplices.push_back( Simplex( {u,v}, d ) );
      }
    }

    return SimplicialComplex( simplices.begin(), simplices.end() );
  };
};

}

}

#endif
