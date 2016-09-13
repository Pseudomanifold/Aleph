#ifndef ALEPH_PERSISTENCE_DIAGRAM_CONVERSION_HH__
#define ALEPH_PERSISTENCE_DIAGRAM_CONVERSION_HH__

#include "PersistenceDiagram.hh"
#include "PersistencePairing.hh"
#include "SimplicialComplex.hh"

#include <map>
#include <vector>

namespace aleph
{

template <
  class Index,
  class Simplex
>
std::vector< PersistenceDiagram<typename Simplex::DataType> > makePersistenceDiagrams( const PersistencePairing<Index>& pairing,
                                                                                       const SimplicialComplex<Simplex>& K )
{
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
    result.push_back( *pair.second );

  return result;
}

}

#endif
