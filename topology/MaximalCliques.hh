#ifndef ALEPH_TOPOLOGY_MAXIMAL_CLIQUES_HH__
#define ALEPH_TOPOLOGY_MAXIMAL_CLIQUES_HH__

#include <iterator>
#include <unordered_set>
#include <vector>

#include "SimplicialComplex.hh"

namespace aleph
{

namespace topology
{

namespace detail
{

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

} // namespace topology

} // namespace aleph

#endif
