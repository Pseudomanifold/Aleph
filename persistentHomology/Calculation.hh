#ifndef ALEPH_PERSISTENT_HOMOLOGY_CALCULATION_HH__
#define ALEPH_PERSISTENT_HOMOLOGY_CALCULATION_HH__

#include "Defaults.hh"

#include "persistenceDiagrams/PersistenceDiagram.hh"
#include "persistenceDiagrams/Calculation.hh"

#include "persistentHomology/PersistencePairing.hh"

#include "SimplicialComplex.hh"
#include "SimplicialComplexConversions.hh"

#include <algorithm>
#include <tuple>
#include <vector>

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

template <
  class ReductionAlgorithm = defaults::ReductionAlgorithm,
  class Representation     = defaults::Representation,
  class Simplex
> std::vector< PersistenceDiagram<typename Simplex::DataType> > calculatePersistenceDiagrams( const SimplicialComplex<Simplex>& K )
{
  auto boundaryMatrix = makeBoundaryMatrix<Representation>( K );
  auto pairing        = calculatePersistencePairing<ReductionAlgorithm>( boundaryMatrix );

  return makePersistenceDiagrams( pairing, K );
}

template <
  class ReductionAlgorithm = defaults::ReductionAlgorithm,
  class Representation     = defaults::Representation,
  class DataType
> PersistenceDiagram<DataType> calculatePersistenceDiagram( const BoundaryMatrix<Representation>& boundaryMatrix,
                                                            const std::vector<DataType>& functionValues )
{
  auto pairing = calculatePersistencePairing<ReductionAlgorithm>( boundaryMatrix );
  return makePersistenceDiagram( pairing, functionValues );
}

}

#endif
