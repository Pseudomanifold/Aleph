#ifndef ALEPH_TOPOLOGY_INTERSECTIONS_HH__
#define ALEPH_TOPOLOGY_INTERSECTIONS_HH__

#include <algorithm>
#include <iterator>
#include <set>
#include <type_traits>
#include <vector>

namespace aleph
{

namespace topology
{

/**
  Intersects two simplices with each other. If the intersection is
  empty, the function returns the empty simplex. The data value of
  the simplex is left unspecified. It should be set by the client.
*/

template <class Simplex> Simplex intersect( const Simplex& s, const Simplex& t ) noexcept
{
  using VertexType = typename Simplex::VertexType;

  std::set<VertexType> sVertices( s.begin(), s.end() );
  std::set<VertexType> tVertices( t.begin(), t.end() );

  std::vector<VertexType> stVertices;
  stVertices.reserve( std::min( s.size(), t.size() ) );

  std::set_intersection( sVertices.begin(), sVertices.end(),
                         tVertices.begin(), tVertices.end(),
                         std::back_inserter( stVertices ) );

  return Simplex( stVertices.begin(), stVertices.end() );
}

/**
  Intersects a simplex with a simplicial complex. To this end, the
  intersection between all simplices in the simplicial complex and
  the current simplex will be calculated.

  The result of this intersection is a set of simplices.
*/

template <class SimplicialComplex, class Simplex> std::set<Simplex> intersect( const SimplicialComplex& K, const Simplex& s ) noexcept
{
  static_assert( std::is_same<Simplex, typename SimplicialComplex::ValueType>::value,
                 "Simplex type and simplicial complex value type must coincide" );

  std::set<Simplex> simplices;

  for( auto&& t : K )
  {
    auto u = intersect(s,t);

    // Only insert non-empty intersections
    if( u )
      simplices.insert( u );
  }

  return simplices;
}

} // namespace topology

} // namespace aleph

#endif
