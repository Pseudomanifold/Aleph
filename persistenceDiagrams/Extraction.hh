#ifndef ALEPH_PERSISTENCE_DIAGRAMS_EXTRACTION_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_EXTRACTION_HH__

#include "persistenceDiagrams/PersistenceDiagram.hh"

#include <algorithm>
#include <vector>

#include <cmath>

namespace aleph
{

/**
  Stores all (signed) persistence values in an output iterator. This is
  a convenience function to simplify common operations. If desired, the
  absolute value of the persistence is being reported.
*/

template <class DataType, class OutputIterator>
void persistence( const PersistenceDiagram<DataType>& D,
                  OutputIterator result,
                  bool useAbsoluteValue = false )
{
  for( auto&& point : D )
    *result++ = useAbsoluteValue ? std::abs( point.persistence() ) : point.persistence();
}

/**
  Similar to aleph::persistence(), but takes the multiplicity of points
  into account. Hence, this function does not offer the option to apply
  the absolute value to the persistence.
*/

template <class DataType, class OutputIterator>
void weightedPersistence( const PersistenceDiagram<DataType>& D,
                          OutputIterator result )
{
  using Point = typename PersistenceDiagram<DataType>::Point;

  std::vector<Point> points;
  points.reserve( D.size() );

  for( auto&& point : D )
    points.push_back( point );

  std::sort( points.begin(), points.end(), [] ( const Point& p, const Point& q )
                                           {
                                             if( p.x() == q.x() )
                                               return p.y() < q.y();
                                             else
                                               return p.x() < q.x();
                                           } );

  // Get all unique points
  std::vector<Point> uniquePoints;
  std::unique_copy( points.begin(), points.end(),
                    std::back_inserter( uniquePoints ) );

  for( auto&& point : uniquePoints )
  {
    auto count  = std::count( points.begin(), points.end(), point );
    auto weight = static_cast<double>( count ) / static_cast<double>( points.size() );

    *result++ = point.persistence() / weight;
  }
}

} // namespace aleph

#endif
