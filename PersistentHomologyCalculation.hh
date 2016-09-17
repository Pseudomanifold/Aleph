#ifndef ALEPH_PERSISTENT_HOMOLOGY_CALCULATION_HH__
#define ALEPH_PERSISTENT_HOMOLOGY_CALCULATION_HH__

#include "Defaults.hh"

#include "PersistenceDiagram.hh"
#include "PersistencePairing.hh"

#include "SimplicialComplex.hh"
#include "SimplicialComplexConversions.hh"

#include "PersistenceDiagramCalculation.hh"
#include "PersistencePairingCalculation.hh"

#include <vector>

namespace aleph
{

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
                                                            const std::vector<double>& functionValues )
{
  auto pairing = calculatePersistencePairing<ReductionAlgorithm>( boundaryMatrix );
  return makePersistenceDiagram( pairing, functionValues );
}

}

#endif
