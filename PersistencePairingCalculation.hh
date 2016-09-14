#ifndef ALEPH_PERSISTENCE_PAIRING_CALCULATION_HH__
#define ALEPH_PERSISTENCE_PAIRING_CALCULATION_HH__

#include "BoundaryMatrix.hh"
#include "PersistencePairing.hh"

#include <algorithm>
#include <tuple>

namespace aleph
{

template <
  class ReductionAlgorithm,
  class Representation
> PersistencePairing<typename Representation::Index> calculatePersistencePairing( const BoundaryMatrix<Representation>& M )
{
  using Index              = typename Representation::Index;
  using PersistencePairing = PersistencePairing<Index>;

  BoundaryMatrix<Representation> B = M;

  ReductionAlgorithm reductionAlgorithm;
  reductionAlgorithm( B );

  PersistencePairing pairing;

  auto numColumns = B.getNumColumns();

  for( Index j = Index(0); j < numColumns; j++ )
  {
    Index i;
    bool valid;

    std::tie( i, valid ) = B.getMaximumIndex( j );
    if( valid )
    {
      auto u = i;
      auto v = j;
      auto w = u;

      if( B.isDualized() )
      {
        u  = numColumns - 1 - v;
        v  = numColumns - 1 - w; // Yes, this is correct!
      }

      pairing.add( u, v );
    }
  }

  std::sort( pairing.begin(), pairing.end() );
  return pairing;
}

}

#endif
