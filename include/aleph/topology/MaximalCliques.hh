#ifndef ALEPH_TOPOLOGY_MAXIMAL_CLIQUES_HH__
#define ALEPH_TOPOLOGY_MAXIMAL_CLIQUES_HH__

#include <iterator>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "math/SparseMatrix.hh"

#include "utilities/UnorderedSetOperations.hh"

#include "SimplicialComplex.hh"


namespace aleph
{

namespace topology
{

namespace detail
{

template <class Simplex> auto adjacencyMatrix( const SimplicialComplex<Simplex>& K ) -> math::SparseBinaryMatrix<typename Simplex::VertexType>
{
  using VertexType = typename Simplex::VertexType;
  using Matrix     = math::SparseBinaryMatrix<VertexType>;

  std::set<VertexType> vertices;
  K.vertices( std::inserter( vertices, vertices.begin() ) );

  // This is required for simplicial complexes in which the vertices do
  // not start at zero.
  std::unordered_map<VertexType, VertexType> vertex_to_id;

  {
    VertexType index = VertexType();
    for( auto&& vertex : vertices )
      vertex_to_id[vertex]  = index++;
  }

  Matrix A( VertexType( vertices.size() ) );

  for( auto itPair = K.range(1); itPair.first != itPair.second; ++itPair.first )
  {
    auto s = *itPair.first;

    // Map vertices (which may be arbitrary natural numbers) to the
    // corresponding column/row in the matrix (which are zero-based
    // and do not have gaps).
    auto u = vertex_to_id.at( s[0] );
    auto v = vertex_to_id.at( s[1] );

    A.set(u,v);
    A.set(v,u);
  }

  A.setIndices( vertices.begin(), vertices.end() );
  return A;
}

template <class VertexType>
void enumerateKoch( std::unordered_set<VertexType>& C,
                    std::unordered_set<VertexType>& I,
                    std::unordered_set<VertexType>& X,
                    std::vector< std::set<VertexType> >& cliques,
                    const math::SparseBinaryMatrix<VertexType>& A )
{
  if( I.empty() && X.empty() )
  {
    std::set<VertexType> newClique;
    for( auto&& c : C )
      newClique.insert( A.getIndex(c) );

    cliques.push_back( newClique );
    return;
  }

  // TODO: Is this the correct way of calculating the set? In the
  // algorithm, it appears that X and I are either both non-empty
  // or empty.
  if( I.empty() )
    return;

  // Pivot selection ---------------------------------------------------

  auto pivot     = *I.begin();
  auto maxDegree = VertexType(0);

  for( auto&& v : I )
  {
    auto degree = VertexType( A.numEntries(v) );
    if( degree > maxDegree )
    {
      pivot     = v;
      maxDegree = degree;
    }
  }

  // Bron--Kerbosch traversal ------------------------------------------

  for( auto it = I.begin(); it != I.end(); )
  {
    auto element = *it;

    // If the currently selected element is a neighbour of the pivot vertex,
    // just move on to another element.
    if( A.get(element, pivot) )
    {
      ++it;
      continue;
    }

    it = I.erase( it );

    std::unordered_set<VertexType> newC = C;
    newC.insert( element );

    std::unordered_set<VertexType> neighbours;
    A.get( element, std::inserter( neighbours, neighbours.begin() ) );

    std::unordered_set<VertexType> newI;
    std::unordered_set<VertexType> newX;

    using namespace aleph::utilities;

    set_intersection( I, neighbours, newI );
    set_intersection( X, neighbours, newX );

    enumerateKoch( newC,
                   newI,
                   newX,
                   cliques,
                   A );

    X.insert( element );
  }
}

template <class VertexType>
void enumerateBronKerbosch( std::unordered_set<VertexType>& C,
                            std::unordered_set<VertexType>& I,
                            std::unordered_set<VertexType>& X,
                            std::vector< std::set<VertexType> >& cliques,
                            const math::SparseBinaryMatrix<VertexType>& A )
{
  if( I.empty() && X.empty() )
  {
    std::set<VertexType> newClique;
    for( auto&& c : C )
      newClique.insert( A.getIndex(c) );

    cliques.push_back( newClique );
    return;
  }

  for( auto it = I.begin(); it != I.end(); )
  {
    auto element = *it;
    it = I.erase( it );

    auto newC = C;
    newC.insert( element );

    std::unordered_set<VertexType> neighbours;
    A.get( element, std::inserter( neighbours, neighbours.begin() ) );

    std::unordered_set<VertexType> newI;
    std::unordered_set<VertexType> newX;

    using namespace aleph::utilities;

    set_intersection( I, neighbours, newI );
    set_intersection( X, neighbours, newX );

    enumerateBronKerbosch( newC,
                           newI,
                           newX,
                           cliques,
                           A );

    X.insert( element );
  }
}

} // namespace detail

/**
  Enumerates all maximal cliques in the given simplicial complex by
  using Koch's modification of the Bron--Kerbosch algorithm for the
  enumeration of cliques.

  Cliques are returned in the form a 2-dimensional vector. For each
  clique, it contains the vertex indices.
*/

template <class Simplex> auto maximalCliquesKoch( const SimplicialComplex<Simplex>& K ) -> std::vector< std::set<typename Simplex::VertexType> >
{
  using VertexType = typename Simplex::VertexType;

  std::unordered_set<VertexType> C;
  std::unordered_set<VertexType> I;
  std::unordered_set<VertexType> X;

  K.vertices( std::inserter( I, I.begin() ) );

  std::vector< std::set<VertexType> > cliques;

  detail::enumerateKoch( C, I, X, cliques, detail::adjacencyMatrix(K) );
  return cliques;
}

template <class Simplex> auto maximalCliquesBronKerbosch( const SimplicialComplex<Simplex>& K ) -> std::vector< std::set<typename Simplex::VertexType> >
{
  using VertexType = typename Simplex::VertexType;

  std::unordered_set<VertexType> C;
  std::unordered_set<VertexType> I;
  std::unordered_set<VertexType> X;

  K.vertices( std::inserter( I, I.begin() ) );

  std::vector< std::set<VertexType> > cliques;

  detail::enumerateBronKerbosch( C, I, X, cliques, detail::adjacencyMatrix(K) );
  return cliques;
}

} // namespace topology

} // namespace aleph

#endif
