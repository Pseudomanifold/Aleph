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

} // namespace filtrations

} // namespace topology

} // namespace aleph

#endif
