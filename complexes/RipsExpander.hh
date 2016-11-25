#ifndef ALEPH_COMPLEXES_RIPS_EXPANDER_HH__
#define ALEPH_COMPLEXES_RIPS_EXPANDER_HH__

#include "Simplex.hh"
#include "SimplicialComplex.hh"

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace aleph
{

namespace complexes
{

template <class Simplex> class RipsExpander
{
public:
  using DataType          = typename Simplex::DataType;
  using VertexType        = typename Simplex::VertexType;

  using SimplicialComplex = SimplicialComplex<Simplex>;
  using SimplexContainer  = std::list<Simplex>;

  // Expansion ---------------------------------------------------------

  SimplicialComplex operator()( const SimplicialComplex& K, unsigned dimension )
  {
    std::set<VertexType> vertices;
    K.vertices( std::inserter( vertices,
                               vertices.begin() ) );

    auto lowerNeighbours = getLowerNeighbours( K );

    std::list<Simplex> simplices;

    for( auto&& vertex : vertices )
    {
      simplices.push_back( Simplex( vertex ) );

      this->addCofaces( simplices.back(),
                        lowerNeighbours.at( vertex ),
                        simplices,
                        dimension );
    }

    return SimplicialComplex( simplices.begin(), simplices.end() );
  }

private:

  using VertexContainer    = std::vector<VertexType>;
  using LowerNeighboursMap = std::unordered_map<VertexType, VertexContainer>;

  static void addCofaces( const Simplex& s,
                          const VertexContainer& neighbours,
                          SimplexContainer& simplices,
                          unsigned dimension )
  {
    if( s.dimension() > dimension )
      return;

    std::unordered_set<VertexType> currentNeighbours( neighbours.begin(),
                                                      neighbours.end() );

    for( auto&& neighbour : neighbours )
    {
      // Create new simplex that contains the new neighbouring vertex as an
      // additional vertex. This increases the dimension by one.
      std::vector<VertexType> vertices( s.begin(), s.end() );
      vertices.push_back( neighbour );

      // TODO:
      //  - Get lower neighbours of current vertex
      //  - Intersect current neighbours with these neighbours
      //  - Add *their* respective cofaces

      simplices.push_back( Simplex( vertices.begin(), vertices.end() ) );
    }
  }

  static LowerNeighboursMap getLowerNeighbours( const SimplicialComplex& K )
  {
    LowerNeighboursMap lowerNeighbours;

    // We only need to traverse the 1-skeleton of the simplicial complex. By
    // adding edges, we automatically fill up all lower neighbours.

    auto&& pair = K.range(1);
    for( auto it = pair.first; it != pair.second; ++it )
    {
      auto&& u = *( it->begin()    ); // first vertex of edge
      auto&& v = *( it->begin() + 1); // second vertex of edge

      if( u < v )
        lowerNeighbours[v].push_back( u );
      else
        lowerNeighbours[u].push_back( v );
    }

    for( auto&& pair : lowerNeighbours )
      std::sort( pair.second.begin(), pair.second.end() );
  }
};

}

}

#endif
