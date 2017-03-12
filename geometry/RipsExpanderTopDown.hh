#ifndef ALEPH_GEOMETRY_RIPS_EXPANDER_TOP_DOWN_HH__
#define ALEPH_GEOMETRY_RIPS_EXPANDER_TOP_DOWN_HH__

#include <list>
#include <vector>

#include "topology/MaximalCliques.hh"

namespace aleph
{

namespace geometry
{

template <class SimplicialComplex> class RipsExpanderTopDown
{
public:
  using Simplex           = typename SimplicialComplex::ValueType;
  using DataType          = typename Simplex::DataType;
  using VertexType        = typename Simplex::VertexType;

  SimplicialComplex operator()( const SimplicialComplex& K, unsigned kMax )
  {
    return this->operator()( K, kMax, 0 );
  }

  SimplicialComplex operator()( const SimplicialComplex& K, unsigned kMax, unsigned kMin )
  {
    (void) kMax;
    (void) kMin;

    auto maximalCliques = aleph::topology::maximalCliquesKoch( K );

    std::list<Simplex> simplices;

    for( auto&& clique : maximalCliques )
      simplices.push_back( Simplex( clique.begin(), clique.end() ) );

    return SimplicialComplex( simplices.begin(), simplices.end() );
  }
};

} // namespace geometry

} // namespace aleph

#endif
