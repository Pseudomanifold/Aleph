#ifndef ALEPH_TOPOLOGY_CLIQUE_GRAPH_HH__
#define ALEPH_TOPOLOGY_CLIQUE_GRAPH_HH__

#include "topology/SimplicialComplex.hh"

#include <map>
#include <stdexcept>
#include <vector>

namespace aleph
{

namespace topology
{

/**
  Given a simplicial complex, extracts its corresponding clique graph. The
  clique graph is defined as the graph in which each node corresponds to a
  k-simplex and an edge connects two nodes if there is a $(k-1)$-face that
  connects the two simplices.

  Note that the graph is represented as a simplicial complex. It makes any
  other operations easier.
*/

template <class Simplex> SimplicialComplex<Simplex> getCliqueGraph( const SimplicialComplex<Simplex>& K, unsigned k )
{
  // Maps k-simplices to their corresponding index in the filtration order of
  // the simplicial complex. This simplifies the creation of edges.
  std::map<Simplex, unsigned> simplexMap;

  // Stores the co-faces of (k-1)-dimensional simplices. This is required for
  // the edge creation. Whenever two (or more) k-simplices appear in this map
  // they will be connected by an edge.
  std::map<Simplex, std::vector<unsigned> > cofaceMap;

  std::vector<Simplex> vertices;
  std::vector<Simplex> edges;

  unsigned index = 0;

  // TODO: Improved traversal possible by taking only the k-dimensional
  // simplices into account here.
  for( auto&& s : K )
  {
    if( s.dimension() == k )
      simplexMap[s] = index;

    ++index;
  }

  for( auto&& pair : simplexMap )
  {
    auto&& simplex = pair.first;
    auto&& index   = pair.second;

    for( auto itFace = simplex.begin_boundary(); itFace != simplex.end_boundary(); ++itFace )
      cofaceMap[ *itFace ].push_back( index );
  }

  // Create vertices ---------------------------------------------------

  for( auto&& pair : simplexMap )
  {
    auto&& simplex = pair.first;
    auto&& index   = pair.second;

    vertices.push_back( Simplex( index, simplex.data() ) );
  }

  // Create edges ------------------------------------------------------

  for( auto&& pair : cofaceMap )
  {
    auto&& indices = pair.second;

    if( indices.size() >= 2 )
    {
      for( std::size_t i = 0; i < indices.size(); i++ )
      {
        auto uIndex = indices[i];

        for( std::size_t j = i+1; j < indices.size(); j++ )
        {
          auto vIndex = indices[j];

          // TODO: What happens if the indices are invalid? Do I need to
          // prepare for this situation explicitly?
          auto&& s = K.at( uIndex );
          auto&& t = K.at( vIndex );

          // TODO: Does it make sense to make this configurable? As of now, I
          // am restricting the weights to the usual maximum filtration, i.e.
          // I am assuming a growth process here.
          auto data = std::max( s.data(), t.data() );

          edges.push_back( Simplex( {uIndex, vIndex}, data ) );
        }
      }
    }
  }

  std::vector<Simplex> simplices;
  simplices.reserve( vertices.size() + edges.size() );
  simplices.insert( simplices.end(), vertices.begin(), vertices.end() );
  simplices.insert( simplices.end(), edges.begin()   , edges.end()    );

  return SimplicialComplex<Simplex>( simplices.begin(), simplices.end() );
}

} // namespace topology

} // namespace aleph

#endif
