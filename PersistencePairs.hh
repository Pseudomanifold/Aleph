#ifndef ALEPH_PERSISTENCE_PAIRS_HH__
#define ALEPH_PERSISTENCE_PAIRS_HH__

#include "BoundaryMatrix.hh"

#include <iostream>
#include <tuple>

namespace aleph
{

template <class ReductionAlgorithm, class Representation> void computePersistencePairs( const BoundaryMatrix<Representation>& M )
{
  using Index = typename Representation::Index;

  BoundaryMatrix<Representation> B = M;

  ReductionAlgorithm reductionAlgorithm;
  reductionAlgorithm( B );

  auto numColumns = B.getNumColumns();

  for( Index j = Index(0); j < numColumns; j++ )
  {
    Index i;
    bool valid;

    std::tie( i, valid ) = B.getMaximumIndex( j );
    if( valid )
      std::cout << "Pair: " << i << "--" << j << std::endl;
  }
}

}

#endif
