#ifndef ALEPH_DISTANCES_HAUSDORFF_HH__
#define ALEPH_DISTANCES_HAUSDORFF_HH__

#include "distances/Infinity.hh"
#include "persistenceDiagrams/PersistenceDiagram.hh"

#include <limits>

namespace aleph
{

namespace distances
{

template <
  class DataType,
  class Distance = InfinityDistance<DataType>
> DataType hausdorffDistance( const PersistenceDiagram<DataType>& D1,
                              const PersistenceDiagram<DataType>& D2,
                              Distance d = Distance() )
{
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

}

}

#endif
