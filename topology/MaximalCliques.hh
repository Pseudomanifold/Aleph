#ifndef ALEPH_TOPOLOGY_MAXIMAL_CLIQUES_HH__
#define ALEPH_TOPOLOGY_MAXIMAL_CLIQUES_HH__

#include <algorithm>
#include <iterator>
#include <set>
#include <unordered_set>
#include <vector>

#include "math/SparseMatrix.hh"

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

  Matrix A( vertices.size() );

  // TODO: This ensures that vertex indices in the simplicial complex
  // start at zero. Otherwise, I would first have to look up an index
  // within the set above.
  for( auto itPair = K.range(1); itPair.first != itPair.second; ++itPair.first )
  {
    auto s = *itPair.first;
    auto u = s[0];
    auto v = s[1];

    A.set(u,v);
    A.set(v,u);
  }

  return A;
}

template <class Simplex>
void enumerateKoch( std::unordered_set<typename Simplex::VertexType>& C,
                    std::unordered_set<typename Simplex::VertexType>& I,
                    std::unordered_set<typename Simplex::VertexType>& X,
                    std::vector< std::vector<typename Simplex::VertexType> >& cliques,
                    const SimplicialComplex<Simplex>& K )
{
  (void) K;

  using VertexType = typename Simplex::VertexType;

  if( I.empty() && X.empty() )
  {
    cliques.push_back( std::vector<VertexType>( C.begin(), C.end() ) );
    return;
  }

  // Pivot selection ---------------------------------------------------

  /*
   * NYI
   */

  // Bron--Kerbosch traversal ------------------------------------------

  /*
   * NYI
   */
}

template <class Simplex>
void enumerateBronKerbosch( std::unordered_set<typename Simplex::VertexType>& C,
                            std::unordered_set<typename Simplex::VertexType>& I,
                            std::unordered_set<typename Simplex::VertexType>& X,
                            std::vector< std::vector<typename Simplex::VertexType> >& cliques,
                            const math::SparseBinaryMatrix<typename Simplex::VertexType>& A )
{
  using VertexType = typename Simplex::VertexType;

  if( I.empty() && X.empty() )
  {
    cliques.push_back( std::vector<VertexType>( C.begin(), C.end() ) );
    return;
  }

  for( auto it = I.begin(); it != I.end(); ++it )
  {
    auto element = *it;
    it = I.erase( it );

    auto newC = C;
    newC.insert( element );

    std::unordered_set<VertexType> neighbours;
    A.get( element, std::inserter( neighbours, neighbours.begin() ) );

    std::unordered_set<VertexType> newI;
    std::unordered_set<VertexType> newX;

    std::set_intersection( I.begin(), I.end(),
                           neighbours.begin(), neighbours.end(),
                           std::inserter( newI,
                                          newI.begin() ) );

    std::set_intersection( X.begin(), X.end(),
                           neighbours.begin(), neighbours.end(),
                           std::inserter( newX,
                                          newX.begin() ) );

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

template <class Simplex> auto maximalCliquesKoch( const SimplicialComplex<Simplex>& K ) -> std::vector< std::vector<typename Simplex::VertexType> >
{
  using VertexType = typename Simplex::VertexType;

  std::unordered_set<VertexType> C;
  std::unordered_set<VertexType> I;
  std::unordered_set<VertexType> X;

  K.vertices( std::inserter( I, I.begin() ) );

  std::vector< std::vector<VertexType> > cliques;

  detail::enumerateKoch( C, I, X, cliques, K );
  return cliques;
}

template <class Simplex> auto maximalCliquesBronKerbosch( const SimplicialComplex<Simplex>& K ) -> std::vector< std::vector<typename Simplex::VertexType> >
{
  using VertexType = typename Simplex::VertexType;

  std::unordered_set<VertexType> C;
  std::unordered_set<VertexType> I;
  std::unordered_set<VertexType> X;

  K.vertices( std::inserter( I, I.begin() ) );

  std::vector< std::vector<VertexType> > cliques;

  detail::enumerateBronKerbosch( C, I, X, cliques, detail::adjacencyMatrix(K) );
  return cliques;
}

} // namespace topology

} // namespace aleph

#endif
