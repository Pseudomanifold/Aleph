#ifndef ALEPH_DISTANCES_NEAREST_NEIGHBOUR_HH__
#define ALEPH_DISTANCES_NEAREST_NEIGHBOUR_HH__

#include <aleph/distances/Infinity.h>>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <algorithm>
#include <limits>
#include <list>
#include <numeric>

namespace aleph
{

namespace distances
{

template <
  class DataType,
  class Distance = InfinityDistance<DataType>
> DataType nearestNeighbourDistance( const PersistenceDiagram<DataType>& D1,
                                     const PersistenceDiagram<DataType>& D2,
                                     Distance d = Distance() )
{
  auto&& oneSidedDistances = [&] ( const PersistenceDiagram<DataType>& D1,
                                   const PersistenceDiagram<DataType>& D2,
                                   Distance d )
  {
    std::list<DataType> distances;

    for( auto&& p1 : D1 )
    {
      DataType infimum = std::numeric_limits<DataType>::max();
      for( auto&& p2 : D2 )
        infimum = std::min( infimum, d( p1, p2 ) );

      distances.push_back( infimum );
    }

    return distances;
  };

  auto&& distances1 = oneSidedDistances( D1, D2, d );
  auto&& distances2 = oneSidedDistances( D2, D1, d );

  return std::accumulate( distances1.begin(), distances1.end(), DataType(0) )
       + std::accumulate( distances2.begin(), distances2.end(), DataType(0) );
}

}

}

#endif
