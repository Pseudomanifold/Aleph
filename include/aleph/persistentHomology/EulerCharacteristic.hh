#ifndef ALEPH_EULER_CHARACTERISTIC_HH__
#define ALEPH_EULER_CHARACTERISTIC_HH__

#include <iterator>
#include <map>

namespace aleph
{

/**
  Calculates the Euler characteristic, i.e. the alternating sum of
  simplex cardinalities, of a given simplicial complex.
*/

template <class SimplicialComplex> long eulerCharacteristic( const SimplicialComplex& K )
{
  // Empty complexes could also be given an invalid characteristic here,
  // but this cannot be expressed through an integer.
  if( K.empty() )
    return 0;

  long chi = 0;

  std::map<std::size_t, std::size_t> cardinality;

  for( auto&& s : K )
  {
    auto d         = s.dimension();
    cardinality[d] = cardinality[d] + 1;
  }

  short s = 1;
  for( std::size_t d = 0; d < K.dimension(); d++ )
  {
    chi += s * cardinality[d];
    s    = s * (-1);
  }

  return chi;
}

/**
  Calculates the Euler characteristic, i.e the alternating sum of Betti
  numbers, from a set of persistence diagrams.
*/

template <class InputIterator> long eulerCharacteristic( InputIterator begin,
                                                         InputIterator end )
{
  if( begin == end )
    return 0;

  using PersistenceDiagram = typename std::iterator_traits<InputIterator>::value_type;

  std::map<std::size_t, std::size_t> betti;

  for( auto it = begin; it != end; ++it )
    betti[ it->dimension() ] = it->betti();

  // Obtain the largest dimension stored by the sequence of persistence
  // diagrams.
  auto D = betti.rbegin()->first;

  long chi = 0;
  short s  = 1;

  for( std::size_t d = 0; d < D; d++ )
  {
    chi += s * betti[d];
    s   *= (-1);
  }

  return chi;
}

/**
  Calculates the persistent Euler characteristic of a sequence of
  persistence diagrams.
*/

template <class InputIterator> auto persistentEulerCharacteristic( InputIterator begin,
                                                                   InputIterator end )
  -> typename std::iterator_traits<InputIterator>::value_type::DataType
{
  using PersistenceDiagram = typename std::iterator_traits<InputIterator>::value_type;
  using DataType           = typename PersistenceDiagram::DataType;

  if( begin == end )
    return DataType();

  DataType chi = DataType();

  for( auto it = begin; it != end; ++it )
  {
    auto&& diagram = *it;
    auto s         = std::pow( -1, it->dimension);

    for( auto&& point : diagram )
      chi += s * point.persistence();
  }

  return chi;

}

} // namespace aleph

#endif
