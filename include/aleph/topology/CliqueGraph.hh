#ifndef ALEPH_TOPOLOGY_CLIQUE_GRAPH_HH__
#define ALEPH_TOPOLOGY_CLIQUE_GRAPH_HH__

#include <aleph/topology/SimplicialComplex.hh>

#include <iterator>
#include <list>
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
  \f$k\f$-simplex and an edge connects two nodes whenever there exists one
  \f$(k-1)\f$-face that connects the two simplices. Edges in the graph are
  weighted, using the *maximum* of the simplex weights of their endpoints;
  this behaviour can be changed by calling an overload of this function.

  @param K Simplicial complex
  @param k Degree of cliques to extract

  Note that the graph is represented as a simplicial complex. It makes any
  further operations easier.
*/

template <class Simplex> SimplicialComplex<Simplex> getCliqueGraph( const SimplicialComplex<Simplex>& K, unsigned k )
{
  using DataType = typename Simplex::DataType;

  return getCliqueGraph( K, k, [] ( DataType a, DataType b ) { return std::max(a,b); } );
}

/**
  Extracts the clique graph using a predefined functor for assigning the
  weights of simplices. The functor requires the following interface:

  \code{.cpp}
  // DataType refers to the data type stored in the simplicial complex,
  // for example `double`.
  DataType Functor::operator()( DataType a, DataType b )
  {
  }
  \endcode

  @param K       Simplicial complex
  @param k       Degree of cliques to extract
  @param functor Functor for assigning weights
*/

template <class Simplex, class Functor> SimplicialComplex<Simplex> getCliqueGraph( const SimplicialComplex<Simplex>& K, unsigned k, Functor functor )
{
  // Stores the co-faces of (k-1)-dimensional simplices. This is required for
  // the edge creation. Whenever two (or more) k-simplices appear in this map
  // they will be connected by an edge.
  std::map<Simplex, std::vector<std::size_t> > cofaceMap;

  std::vector<Simplex> vertices;
  using VertexType = typename Simplex::VertexType;

  {
    // Maps k-simplices to their corresponding index in the filtration order of
    // the simplicial complex. This simplifies the creation of edges.
    std::map<Simplex, std::size_t> simplexMap;

    for( auto itPair = K.range(k); itPair.first != itPair.second; ++itPair.first )
     simplexMap[ *itPair.first ] = K.index( *itPair.first );

    for( auto&& pair : simplexMap )
    {
      auto&& simplex = pair.first;
      auto&& index   = pair.second;

      for( auto itFace = simplex.begin_boundary(); itFace != simplex.end_boundary(); ++itFace )
        cofaceMap[ *itFace ].push_back( index );
    }

    // Create vertices -------------------------------------------------

    vertices.reserve( simplexMap.size() );

    for( auto&& pair : simplexMap )
    {
      auto&& simplex = pair.first;
      auto&& index   = pair.second;

      vertices.push_back( Simplex( VertexType(index), simplex.data() ) );
    }
  }

  // Create edges ------------------------------------------------------

  // Since the number of edges is unknown beforehand, it makes sense to use a
  // list rather than a vector here. Else, the allocation of larger swaths of
  // memory will be problematic.
  std::list<Simplex> edges;

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

          auto&& s    = K.at( uIndex );
          auto&& t    = K.at( vIndex );
          auto data   = functor( s.data(), t.data() );

          edges.push_back( Simplex( {VertexType(uIndex), VertexType(vIndex)}, data ) );
        }
      }
    }
  }

  SimplicialComplex<Simplex> L;
  L.insert( std::make_move_iterator( vertices.begin() ), std::make_move_iterator( vertices.end() ) );
  L.insert( std::make_move_iterator( edges.begin() ), std::make_move_iterator( edges.end() ) );

  return L;
}

} // namespace topology

} // namespace aleph

#endif
