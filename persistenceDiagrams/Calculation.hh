#ifndef ALEPH_PERSISTENCE_DIAGRAM_CALCULATION_HH__
#define ALEPH_PERSISTENCE_DIAGRAM_CALCULATION_HH__

#include "persistenceDiagrams/PersistenceDiagram.hh"
#include "persistentHomology/PersistencePairing.hh"

#include <algorithm>
#include <map>
#include <vector>

namespace aleph
{

/**
  Calculates a set of persistence diagrams from a persistence pairing
  and an associated simplicial complex. The simplicial complex serves
  as a container for looking up the weights for the persistence pairs
  that are stored in the persistence diagram.
*/

template <
  class Index,
  class SimplicialComplex
>
std::vector< PersistenceDiagram<typename SimplicialComplex::ValueType::DataType> > makePersistenceDiagrams( const PersistencePairing<Index>& pairing,
                                                                                                            const SimplicialComplex& K )
{
  using Simplex            = typename SimplicialComplex::ValueType;
  using PersistenceDiagram = PersistenceDiagram<typename Simplex::DataType>;

  std::map<std::size_t, PersistenceDiagram> persistenceDiagrams;

  for( auto&& pair : pairing )
  {
    auto&& i = pair.first;    // Index of creator simplex (always valid)
    auto&& j = pair.second;   // Index of destroyer simplex (may be invalid)

    auto&& s = K.at(i);       // Creator simplex
    auto&& d = s.dimension(); // Current dimension

    if( j < K.size() )
    {
      auto&& t = K.at(j);

      persistenceDiagrams[d].add(
        s.data(), t.data()
      );
    }
    else
      persistenceDiagrams[d].add( s.data() );
  }

  std::vector<PersistenceDiagram> result;
  result.reserve( persistenceDiagrams.size() );

  for( auto&& pair : persistenceDiagrams )
  {
    auto&& diagram = pair.second;
    diagram.setDimension( pair.first );

    result.push_back( diagram );
  }

  std::sort( result.begin(), result.end(),
             [] ( const PersistenceDiagram& D1, const PersistenceDiagram& D2 )
             {
               return D1.dimension() < D2.dimension();
             } );

  return result;
}

/**
  Calculates a persistence diagram from a persistence pairing of a 1D
  function without requiring a representation as a simplicial complex
  for looking up function values.
*/

template <
  class Index,
  class DataType
>
PersistenceDiagram<DataType> makePersistenceDiagram( const PersistencePairing<Index>& pairing,
                                                     const std::vector<DataType>& functionValues )
{
  using PersistenceDiagram = PersistenceDiagram<DataType>;

  PersistenceDiagram D;
  D.setDimension( 0 );

  for( auto&& pair : pairing )
  {
    auto&& i = pair.first;
    auto&& j = pair.second;

    D.add( functionValues.at(i), functionValues.at(j) );
  }

  return D;
}

} // namespace aleph

#endif
