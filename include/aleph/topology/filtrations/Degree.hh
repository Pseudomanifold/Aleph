#ifndef ALEPH_TOPOLOGY_FILTRATIONS_DEGREE_HH__
#define ALEPH_TOPOLOGY_FILTRATIONS_DEGREE_HH__

#include <iterator>
#include <map>
#include <set>

namespace aleph
{

namespace topology
{

namespace filtrations
{

/**
  Calculates *all* vertex degrees of the given simplicial complex. The
  degree of a vertex in a simplicial complex is just the number of its
  co-faces. If the simplicial complex is one-dimensional, degrees in a
  graph-theoretical sense are being calculated.

  @param K      Simplicial complex

  @param result Output iterator for storing the results. The order in
                which degrees are being reported follows the order of
                vertices when calculating them via `vertices()`. Note
                that all degrees are reported as unsigned values.

  @tparam SimplicialComplex Simplicial complex class
  @tparam OutputIterator    Output iterator for storing the results
*/

template <class SimplicialComplex, class OutputIterator> void degrees( const SimplicialComplex& K, OutputIterator result )
{
  using Simplex    = typename SimplicialComplex::ValueType;
  using VertexType = typename Simplex::VertexType;
  using DegreeType = unsigned;

  std::map<VertexType, DegreeType> degrees;

  for( auto&& simplex : K )
  {
    if( simplex.dimension() != 0 )
    {
      for( auto&& vertex : simplex )
        degrees[vertex] += 1;
    }
  }

  for( auto&& pair : degrees )
    *result++ = pair.second;
}

/**
  Calculates the \f$n\f$-degree of a simplex, where \f$n\f$ is the
  dimension of the given simplicial complex. Given another simplex
  \f$\sigma\f$ the \f$n\f$-degree of \f$\sigma\f$ is the number of
  \f$n\f$-simplices that have \f$\sigma\f$ as a face.

  Please refer to the publication *Positively Curved Combinatorial
  3-Manifolds* by Aaron Trout for more details.

  @see http://www.combinatorics.org/ojs/index.php/eljc/article/view/v17i1r49

  @param K Simplicial complex
  @param s Simplex

  @return The \f$n\f$-degree of \p s
*/

template <class SimplicialComplex> unsigned n_degree( const SimplicialComplex& K, const typename SimplicialComplex::ValueType& s )
{
  auto n         = K.dimension();
  auto iterators = K.range( n );

  unsigned d = 0;

  using Simplex    = typename SimplicialComplex::ValueType;
  using VertexType = typename Simplex::VertexType;

  std::set<VertexType> sVertices( s.begin(), s.end() );

  for( auto&& it = iterators.first; it != iterators.second; ++it )
  {
    std::set<VertexType> tVertices( it->begin(), it->end() );
    std::set<VertexType> stVertices;

    std::set_intersection( sVertices.begin(), sVertices.end(),
                           tVertices.begin(), tVertices.end(),
                           std::inserter( stVertices, stVertices.begin() ) );

    if( stVertices.size() == sVertices.size() )
      ++d;
  }

  return d;
}

} // namespace filtrations

} // namespace topology

} // namespace aleph

#endif
