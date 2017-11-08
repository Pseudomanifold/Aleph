#ifndef ALEPH_TOPOLOGY_SPINE_HH__
#define ALEPH_TOPOLOGY_SPINE_HH__

#include <aleph/topology/Intersections.hh>

#include <vector>

namespace aleph
{

namespace topology
{

/**
  Checks whether a simplex in a simplicial complex is principal, i.e.
  whether it is not a proper face of any other simplex in K.
*/

template <class SimplicialComplex, class Simplex> bool isPrincipal( const Simplex& s, const SimplicialComplex& K )
{
  bool principal = true;

  for( auto&& t : K )
  {
    // This check assumes that the simplicial complex is valid, so it
    // suffices to search faces in one dimension _below_ s..
    if( s.dimension()+1 == t.dimension() )
    {
      if( intersect(s,t) == s )
        principal = false;
    }
  }

  return principal;
}

/**
  Checks whether a simplex in a simplicial complex is admissible, i.e.
  the simplex is *principal* and has at least one free face.
*/

template <class SimplicialComplex, class Simplex> bool isAdmissible( const Simplex& s, const SimplicialComplex& K )
{
  if( !isPrincipal(s,K) )
    return false;

  // Check whether a free face exists ----------------------------------

  std::vector<Simplex> faces( s.begin_boundary(), s.end_boundary() );
  std::vector<bool> admissible( faces.size(), false );

  std::size_t i = 0;

  for( auto&& face : faces )
  {
    for( auto&& t : K )
    {
      // This check assumes that the simplicial complex is valid, so it
      // suffices to search faces in one dimension _below_ t.
      if( face.dimension()+1 == t.dimension() && t != s )
      {
        if( intersect(face,t) == face )
        {
          admissible[i] = false;
          break;
        }
      }
    }

    ++i;
  }

  return admissible;
}

/**
  Performs one step of an elementary simplicial collapse in a given
  simplicial complex. The function assumes that the given simplices
  are *valid* for performing the collapse.
*/

template <class SimplicialComplex, class Simplex> SimplicialComplex elementarySimplicialCollapse(
  const Simplex& sigma,
  const Simplex& delta,
  const SimplicialComplex& K )
{
  auto L = K;

  L.remove_without_validation( sigma );
  L.remove_without_validation( delta );
}

} // namespace topology

} // namespace aleph

#endif
