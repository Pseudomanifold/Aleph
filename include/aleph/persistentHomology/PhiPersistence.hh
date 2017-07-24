#ifndef ALEPH_PERSISTENT_HOMOLOGY_PHI_PERSISTENCE_HH__
#define ALEPH_PERSISTENT_HOMOLOGY_PHI_PERSISTENCE_HH__


#include <aleph/topology/Conversions.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <utility>

namespace aleph
{

/**
  Partitions a simplicial complex according to its $\phi$-persistence
  values. This follows the persistent intersection homology algorithm
  in:

    Persistent Intersection Homology
    Paul Bendich and John Harer

  The function expects a simplicial complex $K$ and a function $\phi$
  that determines whether a simplex is proper or not. The function is
  going to create a new simplicial complex. This complex contains all
  proper simplices (in their original order) followed by all improper
  ones.
*/

template <class Simplex, class Function> std::pair<topology::SimplicialComplex<Simplex>, std::size_t> partition( const topology::SimplicialComplex<Simplex>& K, Function phi )
{
  topology::SimplicialComplex<Simplex> L;

  for( auto&& simplex : K )
  {
    if( phi(simplex) )
      L.push_back( simplex );
  }

  auto s = L.size();

  for( auto&& simplex : K )
  {
    if( !phi(simplex) )
      L.push_back( simplex );
  }

  return std::make_pair( L, s );
}

} // namespace aleph

#endif
