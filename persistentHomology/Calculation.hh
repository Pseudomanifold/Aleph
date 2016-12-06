#ifndef ALEPH_PERSISTENT_HOMOLOGY_CALCULATION_HH__
#define ALEPH_PERSISTENT_HOMOLOGY_CALCULATION_HH__

#include "config/Defaults.hh"

#include "persistenceDiagrams/PersistenceDiagram.hh"
#include "persistenceDiagrams/Calculation.hh"

#include "persistentHomology/PersistencePairing.hh"

#include "topology/Conversions.hh"
#include "topology/SimplicialComplex.hh"

#include <algorithm>
#include <tuple>
#include <unordered_set>
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

  std::unordered_set<Index> creators;

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

      // Column j is non-zero. It destroys the feature created by its
      // lowest 1. Hence, i does not remain a creator.
      creators.erase( i );

      if( B.isDualized() )
      {
        u  = numColumns - 1 - v;
        v  = numColumns - 1 - w; // Yes, this is correct!
      }

      pairing.add( u, v );
    }

    // An invalid maximum index indicates that the corresponding column
    // is empty. Hence, we need to think about whether it signifies one
    // feature with infinite persistence.
    else
    {
      // Only add creators that do not belong to the largest dimension
      // of the boundary matrix. Else, there will be a lot of spurious
      // features that cannot be destroyed due to their dimensions.
      if( B.getDimension(j) != B.getDimension() )
        creators.insert( j );
    }
  }

  for( auto&& creator : creators )
  {
    if( B.isDualized() )
      pairing.add( numColumns - 1 - creator );
    else
      pairing.add( creator );
  }

  std::sort( pairing.begin(), pairing.end() );
  return pairing;
}

template <
  class ReductionAlgorithm = defaults::ReductionAlgorithm,
  class Representation     = defaults::Representation,
  class Simplex
> std::vector< PersistenceDiagram<typename Simplex::DataType> > calculatePersistenceDiagrams( const topology::SimplicialComplex<Simplex>& K )
{
  using namespace topology;

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
