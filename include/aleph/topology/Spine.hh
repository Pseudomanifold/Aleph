#ifndef ALEPH_TOPOLOGY_SPINE_HH__
#define ALEPH_TOPOLOGY_SPINE_HH__

#include <aleph/topology/Intersections.hh>

#include <algorithm>
#include <unordered_map>
#include <set>
#include <vector>

namespace aleph
{

namespace topology
{

/**
  Stores coface relationships in a simplicial complex. Given a simplex
  \f$\sigma\f$, the map contains all of its cofaces. Note that the map
  will be updated upon every elementary collapse.
*/

template <class Simplex> using CofaceMap = std::unordered_map<Simplex, std::unordered_set<Simplex> >;

template <class SimplicialComplex> CofaceMap<typename SimplicialComplex::ValueType> buildCofaceMap( const SimplicialComplex& K )
{
  using Simplex   = typename SimplicialComplex::ValueType;
  using CofaceMap = CofaceMap<Simplex>;

  CofaceMap cofaces;
  for( auto&& s : K )
  {
    // Adding an *empty* list of cofaces (so far) for this simplex
    // simplifies the rest of the code because there is no need to
    // check for the existence of a simplex.
    cofaces[s] = {};

    for( auto itFace = s.begin_boundary(); itFace != s.end_boundary(); ++itFace )
      cofaces[ *itFace ].insert( s );
  }

  return cofaces;
}

template <class Simplex> bool isPrincipal( const CofaceMap<Simplex>& cofaces, const Simplex& s )
{
  return cofaces.at( s ).empty();
}

template <class Simplex> Simplex getFreeFace( const CofaceMap<Simplex>& cofaces, const Simplex& s )
{
  if( !isPrincipal( cofaces, s ) )
    return Simplex();

  // Check whether a free face exists ----------------------------------

  for( auto itFace = s.begin_boundary(); itFace != s.end_boundary(); ++itFace )
  {
    auto&& allCofaces = cofaces.at( *itFace );
    if( allCofaces.size() == 1 && allCofaces.find( s ) != allCofaces.end() )
      return *itFace;
  }

  return Simplex();
}

/**
  Checks whether a simplex in a simplicial complex is principal, i.e.
  whether it is not a proper face of any other simplex in K.
*/

template <class SimplicialComplex, class Simplex> bool isPrincipal( const Simplex& s, const SimplicialComplex& K )
{
  // Individual vertices cannot be considered to be principal because
  // they do not have a free face.
  if( s.dimension() == 0 )
    return false;

  bool principal = true;
  auto itPair    = K.range( s.dimension() + 1 );

  for( auto it = itPair.first; it != itPair.second; ++it )
  {
    auto&& t = *it;

    // This check assumes that the simplicial complex is valid, so it
    // suffices to search faces in one dimension _below_ s. Note that
    // the check only has to evaluate the *size* of the intersection,
    // as this is sufficient to determine whether a simplex is a face
    // of another simplex.
    if( sizeOfIntersection(s,t) == s.size() )
      principal = false;
  }

  return principal;
}

/**
  Checks whether a simplex in a simplicial complex is admissible, i.e.
  the simplex is *principal* and has at least one free face.
*/

template <class SimplicialComplex, class Simplex> Simplex isAdmissible( const Simplex& s, const SimplicialComplex& K )
{
  if( !isPrincipal(s,K) )
    return Simplex();

  // Check whether a free face exists ----------------------------------

  std::vector<Simplex> faces( s.begin_boundary(), s.end_boundary() );
  std::vector<bool> admissible( faces.size(), true );

  std::size_t i = 0;
  auto itPair   = K.range( s.dimension() ); // valid range for searches, viz. *all*
                                            // faces in "one dimension up"

  for( auto&& face : faces )
  {
    for( auto it = itPair.first; it != itPair.second; ++it )
    {
      auto&& t = *it;

      // We do not have to check for intersections with the original
      // simplex from which we started---we already know that we are
      // a face.
      if( t != s )
      {
        if( sizeOfIntersection(face,t) == face.size() )
        {
          admissible[i] = false;
          break;
        }
      }
    }

    ++i;
  }

  auto pos = std::find( admissible.begin(), admissible.end(), true );
  if( pos == admissible.end() )
    return Simplex();
  else
    return faces.at( std::distance( admissible.begin(), pos ) );
}

/**
  Calculates all principal faces of a given simplicial complex and
  returns them.
*/

template <class SimplicialComplex> std::unordered_map<typename SimplicialComplex::value_type, typename SimplicialComplex::value_type> principalFaces( const SimplicialComplex& K )
{
  using Simplex = typename SimplicialComplex::value_type;
  auto L        = K;

  std::unordered_map<Simplex, Simplex> admissible;

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

    // FIXME:
    //
    // In case of equal data values, the assignment from above does
    // *not* work and will result in incorrect candidates.
    M = L;

    bool hasFreeFace = false;
    Simplex freeFace = Simplex();

    for( auto itFace = it->begin_boundary(); itFace != it->end_boundary(); ++itFace )
    {
      bool isFace = false;
      for( auto&& simplex : M )
      {
        if( itFace->dimension() + 1 == simplex.dimension() && simplex != *it )
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
      {
        freeFace = *itFace;
        break;
      }
    }

    if( hasFreeFace )
      admissible.insert( std::make_pair( *it, freeFace ) );
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

  return admissible;
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

  // Step 1: obtain initial set of principal faces to start the process
  // of collapsing the complex.
  auto admissible = principalFaces( L );

  // Step 2: collapse until no admissible simplices are left -----------

  while( !admissible.empty() )
  {
    auto s = admissible.begin()->first;
    auto t = admissible.begin()->second;

    L.remove_without_validation( s );
    L.remove_without_validation( t );

    admissible.erase( s );

    // New simplices ---------------------------------------------------
    //
    // Add new admissible simplices that may potentially have been
    // spawned by the removal of s.

    // 1. Add all faces of the principal simplex, as they may
    //    potentially become admissible again.
    std::vector<Simplex> faces( s.begin_boundary(), s.end_boundary() );

    std::for_each( faces.begin(), faces.end(),
      [&t, &L, &admissible] ( const Simplex& s )
      {
        // TODO: rename function
        auto face = isAdmissible( s, L );

        if( t != s && face )
          admissible.insert( std::make_pair( s, face ) );
      }
    );

    // 2. Add all faces othe free face, as they may now themselves
    //    become admissible.
    faces.assign( t.begin_boundary(), t.end_boundary() );

    std::for_each( faces.begin(), faces.end(),
      [&L, &admissible] ( const Simplex& s )
      {
        // TODO: rename function
        auto face = isAdmissible( s, L );

        if( face )
          admissible.insert( std::make_pair( s, face ) );
      }
    );

#if 0
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

        // 1. Add all faces of the principal simplex, as they may
        //    potentially become admissible again.
        std::vector<Simplex> faces( s.begin_boundary(), s.end_boundary() );

        std::for_each( faces.begin(), faces.end(),
          [&t, &L, &admissible] ( const Simplex& s )
          {
            if( t != s && isAdmissible( s, L ) )
              admissible.insert( s );
          }
        );

        // 2. Add all faces othe free face, as they may now themselves
        //    become admissible.
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
#endif

    // The heuristic above is incapable of detecting *all* principal
    // faces of the complex because this may involve searching *all*
    // co-faces. Instead, it is easier to fill up the admissible set
    // here.
    if( admissible.empty() )
      admissible = principalFaces( L );
  }

  return L;
}

} // namespace topology

} // namespace aleph

#endif
