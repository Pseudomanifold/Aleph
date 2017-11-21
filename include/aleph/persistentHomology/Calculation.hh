#ifndef ALEPH_PERSISTENT_HOMOLOGY_CALCULATION_HH__
#define ALEPH_PERSISTENT_HOMOLOGY_CALCULATION_HH__

#include <aleph/config/Defaults.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/Calculation.hh>

#include <aleph/persistentHomology/PersistencePairing.hh>

#include <aleph/topology/Conversions.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <algorithm>
#include <limits>
#include <tuple>
#include <unordered_set>
#include <vector>

namespace aleph
{

/**
  Given a boundary matrix, reduces it and reads off the resulting
  persistence pairing. An optional parameter can be used to force
  the algorithm to stop processing a part of the pairing. This is
  especially relevant for intersection homology, which sets upper
  limits for the validity of an index in the matrix.

  @param M                          Boundary matrix to reduce

  @param includeAllUnpairedCreators Flag indicating whether all unpaired creators should
                                    be included (regardless of their dimension). If set,
                                    this increases the size of the resulting pairing, as
                                    the highest-dimensional columns of the matrix cannot
                                    be reduced any more. The flag is useful, however, in
                                    case one wants to calculate ordinary homology, where
                                    high-dimensional simplices are used for Betti number
                                    calculations.

  @param max                        Optional maximum index after which simplices are not
                                    considered any more. If the pairing of a simplex has
                                    an index larger than the maximum one, such simplices
                                    will not be considered in the pairing. All simplices
                                    are used by default.

  @tparam ReductionAlgorithm Specifies a reduction algorithm to use for reducing
                             the input matrix. Aleph provides a default value in
                             order to simplify the usage of this function.

  @tparam Representation     The representation of the boundary matrix, i.e. how
                             columns are stored. This parameter is automatically
                             determined from the input data.
*/

template <
  class ReductionAlgorithm = aleph::defaults::ReductionAlgorithm,
  class Representation = aleph::defaults::Representation
> PersistencePairing<typename Representation::Index> calculatePersistencePairing( const topology::BoundaryMatrix<Representation>& M,
                                                                                  bool includeAllUnpairedCreators    = false,
                                                                                  typename Representation::Index max = std::numeric_limits<typename Representation::Index>::max() )
{
  using namespace topology;

  using Index              = typename Representation::Index;
  using PersistencePairing = PersistencePairing<Index>;

  BoundaryMatrix<Representation> B = M;

  ReductionAlgorithm reductionAlgorithm;
  reductionAlgorithm( B );

  PersistencePairing pairing;         // resulting pairing
  std::unordered_set<Index> creators; // keeps track of (potential) creators

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

      // Column j is non-zero. It destroys the feature created by its
      // lowest 1. Hence, i does not remain a creator.
      creators.erase( i );

      if( B.isDualized() )
      {
        u  = numColumns - 1 - v;
        v  = numColumns - 1 - w; // Yes, this is correct!
      }

      // u is checked here because it contains the correct index of
      // a simplex with respect to its simplicial complex. Even for
      // a dualized matrix, this index is correctly transformed.
      if( max > numColumns || u < max )
        pairing.add( u, v );
    }

    // An invalid maximum index indicates that the corresponding column
    // is empty. Hence, we need to think about whether it signifies one
    // feature with infinite persistence.
    else
    {
      // Only add creators that do not belong to the largest dimension
      // of the boundary matrix. Else, there will be a lot of spurious
      // features that cannot be destroyed due to their dimensions. If
      // the client wants to have them, however, we let them.
      if(    ( !B.isDualized() && B.getDimension(j) != B.getDimension() )
          || (  B.isDualized() && B.getDimension(j) != Index(0) )
          || includeAllUnpairedCreators )
      {
        creators.insert( j );
      }
    }
  }

  for( auto&& creator : creators )
  {
    auto c = B.isDualized() ? numColumns - 1 - creator : creator;

    // Again, check whether the transformed index value needs to be
    // included in the data. We are not interested in keeping track
    // of simplices that are not allowable (with respect to `max`).
    if( max > numColumns || c < max )
      pairing.add( c );
  }

  std::sort( pairing.begin(), pairing.end() );
  return pairing;
}

/**
  Calculates a set of persistence diagrams from a simplicial complex in
  filtration order, while permitting some additional parameters. Notice
  that this is a *convenience* function that performs *all* conversions
  automatically.

  @param K                          Simplicial complex
  @param dualize                    Indicates that boundary matrix dualization is desired
  @param includeAllUnpairedCreators Indicates that *all* unpaired creators detected during a single pass
                                    of the simplicial complex should be included. This is useful when it
                                    is clear that the simplicial complex models a topological object for
                                    which top-level simplices are meaningful. For Vietoris--Rips complex
                                    calculations, this is usually *not* the case.

  @tparam ReductionAlgorithm Algorithm for reducing the boundary matrix
  @tparam Representation     Representation of the boundary matrix
  @tparam Simplex            Simplex data type (usually inferred from the other parameters)
*/

template <
  class ReductionAlgorithm = defaults::ReductionAlgorithm,
  class Representation     = defaults::Representation,
  class Simplex
> std::vector< PersistenceDiagram<typename Simplex::DataType> > calculatePersistenceDiagrams( const topology::SimplicialComplex<Simplex>& K, bool dualize = true, bool includeAllUnpairedCreators = false )
{
  using namespace topology;

  auto boundaryMatrix = makeBoundaryMatrix<Representation>( K );
  auto pairing        = calculatePersistencePairing<ReductionAlgorithm>( dualize ? boundaryMatrix.dualize() : boundaryMatrix, includeAllUnpairedCreators );

  return makePersistenceDiagrams( pairing, K );
}

/**
  Calculates a persistence diagram from a boundary matrix and a set of
  function values. This function is meant to permit quick calculations
  for one-dimensional functions where a matrix and a vector of values,
  denoting the \f$y\f$-values of the function, are sufficient.

  @param boundaryMatrix Boundary matrix to reduce
  @param functionValues Vector of values to assign to the persistence diagram

  @returns Persistence diagram

  @tparam ReductionAlgorithm Algorithm for reducing the boundary matrix
  @tparam Representation     Representation of the boundary matrix
  @tparam DataType           Data type (usually inferred from the other parameters)
*/

template <
  class ReductionAlgorithm = defaults::ReductionAlgorithm,
  class Representation     = defaults::Representation,
  class DataType
> PersistenceDiagram<DataType> calculatePersistenceDiagram( const topology::BoundaryMatrix<Representation>& boundaryMatrix,
                                                            const std::vector<DataType>& functionValues )
{
  auto pairing = calculatePersistencePairing<ReductionAlgorithm>( boundaryMatrix );
  return makePersistenceDiagram( pairing, functionValues );
}

}

#endif
