#ifndef ALEPH_PERSISTENCE_DIAGRAMS_DISTANCES_POINT_SET_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_DISTANCES_POINT_SET_HH__

/**
  @file  PointSet.hh
  @brief Distance measures based on point sets for persistence diagrams

  The purpose of this file is to collect some distance measures that
  treat a persistence diagram as a simple point set. Hence, parts of
  the topological structure will *not* be considered on purpose. The
  functions based in this file are mainly based on one publication:

  > Distance Measures for Point Sets and Their Computation
  > Thomas Eiter and Heikki Mannila
  > Acta Informatica, Volume 34, Issue 2, pp. 109–133

  @see https://doi.org/10.1007/s002360050075
*/

#include <aleph/geometry/distances/Infinity.hh>

#include <aleph/math/KahanSummation.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <algorithm>
#include <numeric>
#include <vector>

namespace aleph
{

namespace distances
{

/**
  Calculates the sum of minimum distances of points from one diagram to
  the other. The metric is symmetrical because it switches the diagrams
  in its calculations, i.e. all distances to \p D1 are considered first
  and followed by all distances to \p D2.

  @param D1 First persistence diagram
  @param D2 Second persistence diagram

  @param d  Local distance metric for calculating the distance between
            points in the diagram. Defaults to the *infinity distance*
            between points.

  @returns Distance value. Some special cases require extra treatment,
           viz. the case of at least one empty diagram. Here, the best
           course of action is to return `NaN` or, if support for this
           value is lacking, the *maximum* number of the data type. In
           case of both diagrams being empty, zero is returned because
           this is the distance from the empty set to the empty set.
*/

template
<
  class DataType,
  class Distance = aleph::geometry::distances::InfinityDistance<DataType>
> DataType sumOfMinimumDistances( const PersistenceDiagram<DataType>& D1,
                                  const PersistenceDiagram<DataType>& D2,
                                  Distance d = Distance() )
{
  if( D1.empty() && D2.empty() )
    return DataType();
  else if( D1.empty() || D2.empty() )
  {
    if( std::numeric_limits<DataType>::has_infinity )
      return std::numeric_limits<DataType>::infinity();
    else
      return std::numeric_limits<DataType>::max();
  }

  // Distances from points in D1 to D2 ---------------------------------

  std::vector<DataType> distances1;
  distances1.reserve( D1.size() );

  for( auto&& x : D1 )
  {
    DataType distance = std::numeric_limits<DataType>::max();

    for( auto&& y : D2 )
      distance = std::min( distance, d(x,y) );

    distances1.emplace_back( distance );
  }

  // Distances from points in D2 to D1 ---------------------------------

  std::vector<DataType> distances2;
  distances2.reserve( D2.size() );

  for( auto&& y : D2 )
  {
    DataType distance = std::numeric_limits<DataType>::max();

    for( auto&& x : D1 )
      distance = std::min( distance, d(x,y) );

    distances2.emplace_back( distance );
  }

  auto distance1 = aleph::math::accumulate_kahan_sorted( distances1.begin(), distances1.end(), DataType() );
  auto distance2 = aleph::math::accumulate_kahan_sorted( distances2.begin(), distances2.end(), DataType() );

  return DataType(0.5) * (distance1 + distance2);
}

} // namespace distances

} // namespace aleph

#endif
