#ifndef ALEPH_PERSISTENCE_DIAGRAMS_NORMS_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_NORMS_HH__

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/math/KahanSummation.hh>

#include <algorithm>
#include <vector>
#include <stdexcept>

#include <cmath>

namespace aleph
{

/**
  Calculates the total persistence of a given persistence diagram. All
  persistence values will be taken to the $k$th power. Kahan summation
  is used to ensure numerical stability.
*/

template <class DataType> double totalPersistence( const PersistenceDiagram<DataType>& D,
                                                   double k = 2.0,
                                                   bool weighted = false )
{
  aleph::math::KahanSummation<double> result = 0.0;

  if( !weighted )
  {
    for( auto&& point : D )
      result += std::pow( static_cast<double>( std::abs( point.persistence() ) ), k );
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

/** Calculates the $p$-norm of a given persistence diagram. */
template <class DataType> double pNorm( const PersistenceDiagram<DataType>& D,
                                        double p = 2.0,
                                        bool weighted = false )
{
  if( p == 0.0 )
    throw std::runtime_error( "Power must be non-zero" );

  return std::pow( totalPersistence( D, p, weighted ), 1.0 / p );
}

/**
  Calculates the infinity norm of a persistence diagram. The infinity
  norm of a persistence diagram is defined as the maximum persistence
  value stored in the diagram. This function uses absolute values for
  calculating persistence because the result needs to be norm.
*/

template <class DataType> DataType infinityNorm( const PersistenceDiagram<DataType>& D )
{
  if( D.empty() )
    return DataType( 0 );

  std::vector<DataType> persistenceValues;

  for( auto&& point : D )
    persistenceValues.push_back( std::abs( point.persistence() ) );

  return *std::max_element( persistenceValues.begin(), persistenceValues.end() );
}

} // namespace aleph

#endif
