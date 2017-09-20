#ifndef ALEPH_PERSISTENCE_DIAGRAMS_DISTANCES_NEAREST_NEIGHBOUR_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_DISTANCES_NEAREST_NEIGHBOUR_HH__

#include <aleph/geometry/distances/Infinity.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/math/KahanSummation.hh>

#include <algorithm>
#include <limits>
#include <list>
#include <numeric>

namespace aleph
{

namespace distances
{

/**
  Calculates a pseudo-distance by assessing the distance of every point
  in the diagram to its nearest neighbour, measured by some distance on
  the persistence diagram. To make this symmetrical, one-sided distance
  calculations are performed for every point and their sum is returned.

  The purpose of this function is to yield suitable *baselines* for the
  actual distance between two persistence diagrams.

  @param D1 First persistence diagram
  @param D2 Second persistence diagram
  @param d  Distance functor; this argument is required to permit type
            detection by the compiler

  @returns Sum of nearest neighbour estimates
*/

template <
  class DataType,
  class Distance = aleph::geometry::distances::InfinityDistance<DataType>
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

  using namespace math;

  auto sum =   accumulate_kahan_sorted( distances1.begin(), distances1.end(), DataType(0) )
             + accumulate_kahan_sorted( distances2.begin(), distances2.end(), DataType(0) );

  return sum / 2;
}

} // namespace distances

} // namespace aleph

#endif
