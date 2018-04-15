#ifndef ALEPH_EULER_CHARACTERISTIC_HH__
#define ALEPH_EULER_CHARACTERISTIC_HH__

#include <map>

namespace aleph
{

/**
  Calculates the Euler characteristic, i.e. the alternating sum of
  simplex cardinalities, of a given simplicial complex.
*/

template <class SimplicialComplex> eulerCharacteristic( const SimplicialComplex& K )
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

} // namespace aleph

#endif
