#ifndef ALEPH_PERSISTENCE_DIAGRAMS_NORMS_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_NORMS_HH__

#include "persistenceDiagrams/PersistenceDiagram.hh"

#include <algorithm>
#include <vector>
#include <stdexcept>

#include <cmath>

namespace aleph
{

template <class DataType> double totalPersistence( const PersistenceDiagram<DataType>& D,
                                                   double k = 2.0 )
{
  double result = 0.0;

  for( auto&& point : D )
    result = result + std::pow( static_cast<double>( point.persistence() ), k );

  return result;
}

template <class DataType> double pNorm( const PersistenceDiagram<DataType>& D,
                                        double p = 2.0 )
{
  if( p == 0.0 )
    throw std::runtime_error( "Power must be non-zero" );

  return std::pow( totalPersistence( D, p ), 1.0 / p );
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

}

#endif
