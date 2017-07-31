#ifndef ALEPH_PERSISTENCE_DIAGRAMS_DISTANCES_HAUSDORFF_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_DISTANCES_HAUSDORFF_HH__

#include <aleph/geometry/distances/Infinity.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <limits>

namespace aleph
{

namespace distances
{

/**
  Calculates the Hausdorff distance between two persistence diagrams,
  i.e. the Hausdorff distance between their corresponding point sets
  treated as 2D sets.

  By default, the infinity distance ($L_\infty$) is used.

  There are two special cases handled by this function:

  - If both persistence diagrams are empty, a distance of zero will be
    returned. This is required in order to be consistent with a metric
    in mathematics.

  - If exactly one persistence diagram is empty, a distance of +inf is
    returned. This indicates a potentially problematic situation. When
    a given data type does not support positive infinity, its positive
    maximum value is returned.
*/

template <
  class DataType,
  class Distance = InfinityDistance<DataType>
> DataType hausdorffDistance( const PersistenceDiagram<DataType>& D1,
                              const PersistenceDiagram<DataType>& D2,
                              Distance d = Distance() )
{
  if( D1.empty() && D2.empty() )
    return DataType();
  else if( D1.empty() ^ D2.empty() )
  {
    if( std::numeric_limits<DataType>::has_infinity )
      return std::numeric_limits<DataType>::infinity();
    else
      return std::numeric_limits<DataType>::max();
  }

  using PersistenceDiagram = PersistenceDiagram<DataType>;
  using Point              = typename PersistenceDiagram::Point;

  auto&& infimumDistance = [&] ( const Point& p,
                                 const PersistenceDiagram& D,
                                 Distance d )
  {
    DataType infimum = std::numeric_limits<DataType>::max();

    for( auto&& q : D )
      infimum = std::min( infimum, d( p, q ) );

    return infimum;
  };

  DataType supremum1 = std::numeric_limits<DataType>::lowest();

  for( auto&& p : D1 )
    supremum1 = std::max( supremum1, infimumDistance( p, D2, d ) );

  DataType supremum2 = std::numeric_limits<DataType>::lowest();

  for( auto&& p : D2 )
    supremum2 = std::max( supremum2, infimumDistance( p, D1, d ) );

  return std::max( supremum1, supremum2 );
}

} // namespace distances

} // namespace aleph

#endif
