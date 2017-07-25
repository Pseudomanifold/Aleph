#ifndef ALEPH_TOPOLOGY_INTERSECTIONS_HH__
#define ALEPH_TOPOLOGY_INTERSECTIONS_HH__

#include <algorithm>
#include <iterator>
#include <set>
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

} // namespace topology

} // namespace aleph

#endif
