#ifndef ALEPH_TOPOLOGY_SPINE_HH__
#define ALEPH_TOPOLOGY_SPINE_HH__

#include <aleph/topology/Intersections.hh>

#include <unordered_set>
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
  std::vector<bool> admissible( faces.size(), true );

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

  return std::find( admissible.begin(), admissible.end(), true ) != admissible.end();
}

/**
  Checks whether a pair of a simplex and its face are admissible, i.e.
  the simplex is *principal* and the face is free.
*/

template <class SimplicialComplex, class Simplex> bool isAdmissible( const Simplex& sigma, const Simplex& delta, const SimplicialComplex& K )
{
  if( !isPrincipal(sigma,K) )
    return false;

  // Check whether the face is free ------------------------------------

  bool admissible = true;

  // TODO: the traversal could be sped up by only taking at a look at
  // the *relevant* dimension of the simplicial complex.
  for( auto&& s : K )
  {
    // This check assumes that the simplicial complex is valid, so it
    // suffices to search faces in one dimension _below_ delta.
    if( delta.dimension()+1 == s.dimension() && s != sigma )
    {
      if( intersect(delta, s) == delta )
      {
        admissible = false;
        break;
      }
    }
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
  auto L        = K;

  std::unordered_set<Simplex> admissible;

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
      bool isFace = false;
      for( auto&& simplex : M )
      {
        if( itFace->dimension() + 1 == simplex.dimension() )
        {
          // The current face must *not* be a face of another simplex in
          // the simplicial complex.
          if( intersect( *itFace, simplex ) == *itFace )
          {
            isFace = true;
            break;
          }
        }
      }

      hasFreeFace = !isFace;
      if( hasFreeFace )
        break;
    }

    if( hasFreeFace )
      admissible.insert( *it );
  }

  // Step 2: determine principality ------------------------------------
  //
  // All simplices that are faces of higher-dimensional simplices are
  // now removed from the map of admissible simplices.

  for( auto&& s : L )
  {
    for( auto itFace = s.begin_boundary(); itFace != s.end_boundary(); ++itFace )
      admissible.erase( *itFace );
  }

  // Step 3: collapse until no admissible simplices are left -----------

  while( !admissible.empty() )
  {
    auto s           = *admissible.begin();
    bool hasFreeFace = false;

    // TODO: this check could be simplified by *storing* the free face
    // along with the given simplex
    for( auto itFace = s.begin_boundary(); itFace != s.end_boundary(); ++itFace )
    {
      auto t = *itFace;

      if( isAdmissible( s, t, L ) )
      {
        L.remove_without_validation( s );
        L.remove_without_validation( t );

        admissible.erase( s );

        // New simplices -----------------------------------------------
        //
        // Add new admissible simplices that may potentially have been
        // spawned by the removal of s.

        std::vector<Simplex> faces( s.begin_boundary(), s.end_boundary() );

        std::for_each( faces.begin(), faces.end(),
          [&t, &L, &admissible] ( const Simplex& s )
          {
            if( t != s && isAdmissible( s, L ) )
              admissible.insert( s );
          }
        );

        faces.assign( t.begin_boundary(), t.end_boundary() );

        std::for_each( faces.begin(), faces.end(),
          [&L, &admissible] ( const Simplex& s )
          {
            if( isAdmissible( s, L ) )
              admissible.insert( s );
          }
        );

        hasFreeFace = true;
        break;
      }
    }

    // The admissible simplex does not have a free face, so it must not
    // be used.
    if( !hasFreeFace )
      admissible.erase( s );
  }

  return L;
}

} // namespace topology

} // namespace aleph

#endif
