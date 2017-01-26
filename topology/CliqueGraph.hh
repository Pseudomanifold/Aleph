#ifndef ALEPH_TOPOLOGY_CLIQUE_GRAPH_HH__
#define ALEPH_TOPOLOGY_CLIQUE_GRAPH_HH__

#include "topology/CliqueGraph.hh"

#include <map>
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

template <class Simplex> getCliqueGraph( const SimplicialComplex& K, unsigned k )
{
  // Maps k-simplices to their corresponding index in the filtration order of
  // the simplicial complex. This simplifies the creation of edges.
  std::map<Simplex, unsigned> simplexMap;

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

  // Create vertices ---------------------------------------------------

  for( auto&& pair : simplexMap )
  {
    auto&& s = pair.first;
    auto&& v = pair.second;

    vertices.push_back( Simplex( v, s.data() ) );
  }

  return SimplicialComplex( vertices.begin(), vertices.end() );
}

} // namespace topology

} // namespace aleph

#endif
