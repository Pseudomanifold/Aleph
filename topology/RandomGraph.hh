#ifndef ALEPH_TOPOLOGY_RANDOM_GRAPH_HH__
#define ALEPH_TOPOLOGY_RANDOM_GRAPH_HH__

#include "Simplex.hh"
#include "SimplicialComplex.hh"

#include <random>
#include <vector>

namespace aleph
{

namespace topology
{

/**
  Generates an Erdős--Rényi graph with n vertices and a link probability
  of p. Note that the graph will be returned as an unweighted simplicial
  complex.
*/

auto generateErdosRenyiGraph( unsigned n, double p ) -> SimplicialComplex< Simplex<short, unsigned> >
{
  using S = Simplex<short, unsigned>;
  using K = SimplicialComplex<S>;

  std::vector<S> simplices;

  std::random_device rd;
  std::mt19937 mt( rd() );
  std::uniform_real_distribution<> distribution( 0.0, 1.0 );

  for( unsigned i = 0; i < n; i++ )
    simplices.push_back( S( i ) );

  for( unsigned u = 0; u < n; u++ )
    for( unsigned v = u+1; v < n; v++ )
      if( distribution( mt ) < p )
        simplices.push_back( S( {u,v} ) );

  return K( simplices.begin(), simplices.end() );
}


} // namespace topology

} // namespace aleph

#endif
