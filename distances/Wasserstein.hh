#ifndef ALEPH_DISTANCES_WASSERSTEIN_HH__
#define ALEPH_DISTANCES_WASSERSTEIN_HH__

#include "Infinity.hh"
#include "PersistenceDiagram.hh"

#include "detail/Munkres.hh"

#include <algorithm>
#include <limits>
#include <stdexcept>

#include <cmath>

namespace aleph
{

namespace distances
{

template <
  class DataType,
  class Distance = InfinityDistance<DataType>
> DataType wassersteinDistance( const PersistenceDiagram<DataType>& D1,
                                const PersistenceDiagram<DataType>& D2,
                                DataType power = DataType( 1 ) )
{
  if( D1.dimension() != D2.dimension() )
    throw std::runtime_error( "Dimensions do not coincide" );

  auto size = D1.size() + D2.size();

  detail::Matrix<DataType> costs( size );

  using IndexType = decltype( costs.n() );

  IndexType row = IndexType();
  IndexType col = IndexType();

  Distance dist;

  // Regular block of matrix -------------------------------------------
  //
  // This block stores the distances between individual points of the
  // respective persistence diagrams.

  for( auto&& p1 : D1 )
  {
    col = IndexType();

    for( auto&& p2 : D2 )
    {
      auto d = std::pow( dist( p1, p2 ), power );

      costs( row, col ) = d;

      // This block needs to be completely zero, as it contains the distances
      // between the individual orthogonal projections of the points.
      costs( col + D1.size(), row + D2.size() ) = DataType();

      ++col;
    }

    ++row;
  }

  // Orthogonal projection of the first diagram ------------------------

  row = IndexType();
  col = D2.size();

  for( auto&& p1 : D1 )
  {
    row = IndexType();

    for( auto&& p2 : D1 )
    {
      DataType d = DataType();

      if( p1 == p2 )
      {
        // TODO: Orthogonal distance
      }
      else
        d = std::numeric_limits<DataType>::max();

      costs( row, col ) = d;

      ++row;
    }

    ++col;
  }

  // Orthogonal projection of the second diagram -----------------------

  row = D1.size();
  col = IndexType();

  for( auto&& p1 : D2 )
  {
    col = IndexType();

    for( auto&& p2 : D2 )
    {
      DataType d = DataType();

      if( p1 == p2 )
      {
        // TODO: Orthogonal distance
      }
      else
        d = std::numeric_limits<DataType>::max();

      costs( row, col ) = d;

      ++col;
    }

    ++row;
  }

  // Assignment problem solving ----------------------------------------

  detail::Munkres<DataType> solver( costs );

  auto M              = solver();
  DataType totalCosts = DataType();

  for( row = IndexType(); row < M.n(); row++ )
  {
    for( col = IndexType(); col < M.n(); col++ )
    {
      if( M( row, col ) == IndexType() )
        totalCosts += costs( row, col );
    }
  }

  return std::pow( totalCosts, 1 / power );
}

}

}

#endif
