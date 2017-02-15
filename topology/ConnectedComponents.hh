#ifndef ALEPH_TOPOLOGY_CONNECTED_COMPONENTS_HH__
#define ALEPH_TOPOLOGY_CONNECTED_COMPONENTS_HH__

#include "topology/UnionFind.hh"

#include <iterator>
#include <set>

namespace aleph
{

namespace topology
{

/**
  Calculates 'ordinary' connected components of a simplicial complex,
  resulting in a 'Union--Find' data structure that contains simplices
  that are part of the same connected components.

  A client should use the \c UnionFind::roots() function in order to
  get all simplices that create a connected component. Subsequently,
  one should use \c UnionFind::get() with each root vertex to obtain
  all the simplices that make up the given connected component.
*/

template <class SimplicialComplex> UnionFind<typename SimplicialComplex::ValueType::VertexType> calculateConnectedComponents( const SimplicialComplex& K )
{
  using Simplex = typename SimplicialComplex::ValueType;
  using Vertex  = typename Simplex::VertexType;

  std::set<Vertex> vertices;
  K.vertices( std::inserter( vertices, vertices.begin() ) );

  UnionFind<Vertex> uf( vertices.begin(), vertices.end() );

  for( auto itPair = K.range(1); itPair.first != itPair.second; ++itPair.first )
  {
    auto s = *itPair.first;
    auto u = s[0];
    auto v = s[1];

    // Follow the 'elder rule' for consistency reasons. It makes no
    // difference for the output, however.
    if( u < v )
      uf.merge( v, u );
    else
      uf.merge( u, v );
  }

  return uf;
}

} // namespace "topology"

} // namespace "aleph"

#endif
