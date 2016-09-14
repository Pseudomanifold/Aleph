#ifndef ALEPH_PERSISTENT_HOMOLOGY_CALCULATION_HH__
#define ALEPH_PERSISTENT_HOMOLOGY_CALCULATION_HH__

#include "Defaults.hh"

#include "PersistenceDiagram.hh"
#include "PersistencePairing.hh"

#include "SimplicialComplex.hh"
#include "SimplicialComplexConversions.hh"

#include "PersistenceDiagramCalculation.hh"
#include "PersistencePairingCalculation.hh"

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

}

#endif
