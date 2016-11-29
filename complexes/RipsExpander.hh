#ifndef ALEPH_COMPLEXES_RIPS_EXPANDER_HH__
#define ALEPH_COMPLEXES_RIPS_EXPANDER_HH__

#include <algorithm>
#include <list>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace aleph
{

namespace complexes
{

template <class SimplicialComplex> class RipsExpander
{
public:
  using Simplex           = typename SimplicialComplex::ValueType;
  using DataType          = typename Simplex::DataType;
  using VertexType        = typename Simplex::VertexType;

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

  // Weight assignment -------------------------------------------------

  SimplicialComplex assignMaximumWeight( const SimplicialComplex& K, unsigned minDimension = 1 )
  {
    SimplicialComplex S;

    for( auto s : K )
    {
      // Re-calculate the weight of the simplex because its
      // dimensionality requirement is not satisfied
      if( s.dimension() > minDimension )
      {
        auto w = s.data();

        for( auto itFace = s.begin_boundary(); itFace != s.end_boundary(); ++itFace )
        {
          auto itFaceInK = K.find( *itFace );
          if( itFaceInK != K.end() )
            w = std::max( w, itFaceInK->data() );
        }

        s.setData( w );
      }

      // TODO: Not sure whether this is the best way of solving it;
      // should I expect a generic simplicial complex to have this
      // function?
      S.push_back_without_validation( s );
    }

    return S;
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
