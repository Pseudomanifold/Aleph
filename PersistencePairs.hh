#ifndef ALEPH_PERSISTENCE_PAIRS_HH__
#define ALEPH_PERSISTENCE_PAIRS_HH__

#include "Dualization.hh"
#include "BoundaryMatrix.hh"
#include "PersistencePairing.hh"

#include <algorithm>
#include <iostream>
#include <tuple>

namespace aleph
{

template <
  class ReductionAlgorithm,
  class Representation
> PersistencePairing<typename Representation::Index> computePersistencePairs( const BoundaryMatrix<Representation>& M, bool dualize = false )
{
  using Index              = typename Representation::Index;
  using PersistencePairing = PersistencePairing<Index>;

  BoundaryMatrix<Representation> B = M;
  if( dualize )
    dualizeTrivial( B );

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

      if( dualize )
      {
        u  = numColumns - 1 - v;
        v  = numColumns - 1 - w; // Yes, this is correct!
      }

      std::cout << "Pair: " << u << "--" << v << std::endl;

      pairing.add( u, v );
    }
  }

  std::sort( pairing.begin(), pairing.end() );
  return pairing;
}

}

#endif
