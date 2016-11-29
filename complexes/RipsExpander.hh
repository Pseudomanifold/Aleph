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

      if( lowerNeighbours.find( vertex ) != lowerNeighbours.end() )
      {
        addCofaces( simplices.back(),
                    lowerNeighbours,
                    lowerNeighbours.at( vertex ),
                    simplices,
                    dimension );
      }
    }

    return SimplicialComplex( simplices.begin(), simplices.end() );
  }

private:

  using VertexContainer    = std::unordered_set<VertexType>;
  using LowerNeighboursMap = std::unordered_map<VertexType, VertexContainer>;

  static void addCofaces( const Simplex& s,
                          const LowerNeighboursMap& lowerNeighboursMap,
                          const VertexContainer& neighbours,
                          SimplexContainer& simplices,
                          unsigned dimension )
  {
    if( s.dimension() > dimension )
      return;

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
      //  - Set weight

      simplices.push_back( Simplex( vertices.begin(), vertices.end() ) );

      if( lowerNeighboursMap.find( neighbour ) != lowerNeighboursMap.end() )
      {
        auto lowerNeighbours  = lowerNeighboursMap.at( neighbour );
        auto commonNeighbours = intersect( lowerNeighbours, neighbours );

        addCofaces( simplices.back(),
                    lowerNeighboursMap,
                    commonNeighbours,
                    simplices,
                    dimension );
      }
    }
  }

  static std::unordered_set<VertexType> intersect( const std::unordered_set<VertexType>& U,
                                                   const std::unordered_set<VertexType>& V )
  {
    auto uPtr = &U;
    auto vPtr = &V;

    std::unordered_set<VertexType> result;

    if( uPtr->size() > vPtr->size() )
      std::swap( uPtr, vPtr );

    for( auto&& u : *uPtr )
      if( vPtr->find( u ) != vPtr->end() )
        result.insert( u );

    return result;
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
        lowerNeighbours[v].insert( u );
      else
        lowerNeighbours[u].insert( v );
    }

    return lowerNeighbours;
  }
};

}

}

#endif
