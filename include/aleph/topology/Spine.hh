#ifndef ALEPH_TOPOLOGY_SPINE_HH__
#define ALEPH_TOPOLOGY_SPINE_HH__

#include <aleph/topology/Intersections.hh>

#include <unordered_map>
#include <vector>

// FIXME: remove after debugging
#include <iostream>

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

/**
  Performs an iterated elementary simplicial collapse until *all* of the
  admissible simplices have been collapsed. This leads to the *spine* of
  the simplicial complex.

  @see S. Matveev, "Algorithmic Topology and Classification of 3-Manifolds"
*/

template <class SimplicialComplex> SimplicialComplex spine( const SimplicialComplex& K )
{
  using Simplex = typename SimplicialComplex::value_type;

  auto L = K;
  L.sort();

  std::unordered_map<Simplex, bool> admissible;

  // Step 1: determine free faces --------------------------------------
  //
  // This first checks which simplices have at least one free face,
  // meaning that they may be potentially admissible.

  for( auto it = L.begin(); it != L.end(); ++it )
  {
    if( it->dimension() == 0 )
      continue;

    // The range of the complex M is sufficient because we have
    // already encountered all lower-dimensional simplices that
    // precede the current one given by `it`.
    //
    // This complex will be used for testing free faces.
    SimplicialComplex M( L.begin(), it );

    bool hasFreeFace = false;
    for( auto itFace = it->begin_boundary(); itFace != it->end_boundary(); ++itFace )
    {
      for( auto&& simplex : M )
      {
        if( itFace->dimension() + 1 == simplex.dimension() )
        {
          // The current face must *not* be a face of another simplex in
          // the simplicial complex.
          if( intersect( *itFace, simplex ) != *itFace )
          {
            hasFreeFace = true;
            break;
          }

          // This simplex can *never* be a free face, so let's precede
          // to the next face.
          else
            break;
        }
      }
    }

    if( hasFreeFace )
      admissible[ *it ] = true;
  }

  // Step 2: determine principality ------------------------------------
  //
  // All simplices that are faces of higher-dimensional simplices are
  // now removed from the map of admissible simplices.

  for( auto&& s : L )
  {
    for( auto itFace = s.begin_boundary(); itFace != s.end_boundary(); ++itFace )
      admissible[ *itFace ] = false;
  }

  std::cerr << "ADMISSIBLE SIMPLICES:\n";

  for( auto&& pair : admissible )
    if( pair.second )
      std::cerr << pair.first << "\n";

  return L;
}

} // namespace topology

} // namespace aleph

#endif
