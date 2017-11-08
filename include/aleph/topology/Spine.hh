#ifndef ALEPH_TOPOLOGY_SPINE_HH__
#define ALEPH_TOPOLOGY_SPINE_HH__

#include <aleph/topology/Intersections.hh>

namespace aleph
{

namespace topology
{

template <class SimplicialComplex, class Simplex> bool isPrincipal( const Simplex& s, const SimplicialComplex& K )
{
  bool principal = true;

  for( auto&& t : K )
  {
    if( s.dimension() < t.dimension() )
    {
      if( intersect(s,t) == s )
        principal = false;
    }
  }

  return principal;

}

template <class SimplicialComplex, class Simplex> bool isAdmissible( const Simplex& s, const SimplicialComplex& K )
{
  bool admissible = true;

  for( auto&& t : K )
  {
  }

  return admissible;
}

} // namespace topology

} // namespace aleph

#endif
