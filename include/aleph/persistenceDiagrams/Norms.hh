#ifndef ALEPH_PERSISTENCE_DIAGRAMS_NORMS_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_NORMS_HH__

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <algorithm>
#include <vector>
#include <stdexcept>

#include <cmath>

namespace aleph
{

template <class DataType> double totalPersistence( const PersistenceDiagram<DataType>& D,
                                                   double k = 2.0,
                                                   bool weighted = false )
{
  double result = 0.0;

  if( !weighted )
  {
    for( auto&& point : D )
      result = result + std::pow( static_cast<double>( point.persistence() ), k );
  }
  else
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

    std::vector<Point> uniquePoints;

    std::unique_copy( points.begin(), points.end(),
                      std::back_inserter( uniquePoints ) );

    std::vector<unsigned> counts;
    counts.reserve( uniquePoints.size() );

    for( auto&& point : uniquePoints )
    {
      auto count = static_cast<unsigned>( std::count( points.begin(), points.end(), point ) );
      counts.push_back( count );
    }

    auto itPoint = uniquePoints.begin();
    auto itCount = counts.begin();

    for( ; itPoint != uniquePoints.end() && itCount != counts.end(); ++itPoint, ++itCount )
    {
      double weight = *itCount / static_cast<double>( points.size() );
      result       += weight * std::pow( static_cast<double>( itPoint->persistence() ), k );
    }
  }

  return result;
}

template <class DataType> double pNorm( const PersistenceDiagram<DataType>& D,
                                        double p = 2.0,
                                        bool weighted = false )
{
  if( p == 0.0 )
    throw std::runtime_error( "Power must be non-zero" );

  return std::pow( totalPersistence( D, p, weighted ), 1.0 / p );
}

template <class DataType> DataType infinityNorm( const PersistenceDiagram<DataType>& D )
{
  if( D.empty() )
    return DataType( 0 );

  std::vector<DataType> persistenceValues;

  for( auto&& point : D )
    persistenceValues.push_back( point.persistence() );

  return *std::max_element( persistenceValues.begin(), persistenceValues.end() );
}

} // namespace aleph

#endif
