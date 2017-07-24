#ifndef ALEPH_PERSISTENT_HOMOLOGY_PHI_PERSISTENCE_HH__
#define ALEPH_PERSISTENT_HOMOLOGY_PHI_PERSISTENCE_HH__

#include <aleph/topology/SimplicialComplex.hh>

namespace aleph
{

/**
  Partitions a simplicial complex according to its $\phi$-persistence
  values. This follows the persistent intersection homology algorithm
  in:

    Persistent Intersection Homology
    Paul Bendich and John Harer

  The function expects a simplicial complex $K$ and a function $\phi$
  that determines whether a simplex is proper or not.

  TODO: finish documentation; what is a good return type?
*/

template <class Simplex, class Function> void partition( const topology::SimplicialComplex<Simplex>& K, Function phi )
{
}

} // namespace aleph

#endif
