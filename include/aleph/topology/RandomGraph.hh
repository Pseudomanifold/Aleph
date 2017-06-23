#ifndef ALEPH_TOPOLOGY_RANDOM_GRAPH_HH__
#define ALEPH_TOPOLOGY_RANDOM_GRAPH_HH__

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <algorithm>
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

/**
  Generates a weighted random graph with n vertices and a link
  probability of p. In contrast to Erdős--Rényi graphs, here a
  weight is assigned according to a number of Bernoulli trials
  with success probability p.
*/

auto generateWeightedRandomGraph( unsigned n, double p ) -> SimplicialComplex< Simplex<unsigned, unsigned> >
{
  using S = Simplex<unsigned, unsigned>;
  using K = SimplicialComplex<S>;

  std::vector<S> simplices;

  std::random_device rd;
  std::mt19937 mt( rd() );
  std::uniform_real_distribution<> distribution( 0.0, 1.0 );

  for( unsigned i = 0; i < n; i++ )
    simplices.push_back( S( i ) );

  for( unsigned u = 0; u < n; u++ )
  {
    for( unsigned v = u+1; v < n; v++ )
    {
      unsigned w = 0;

      // Repeated Bernoulli trials: until the first failure occurs, the
      // edge weight is increased.
      while( distribution( mt ) < p )
        ++w;

      if( w != 0 )
        simplices.push_back( S( {u,v}, w ) );
    }
  }

  std::sort( simplices.begin(), simplices.end(), aleph::topology::filtrations::Data<S>() );

  return K( simplices.begin(), simplices.end() );

}

} // namespace topology

} // namespace aleph

#endif
