#ifndef ALEPH_PERSISTENCE_DIAGRAMS_MEAN_HH__
#define ALEPH_PERSISTENCE_DIAGRAMS_MEAN_HH__

#include <aleph/distances/Infinity.hh>

#include <aleph/math/KahanSummation.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistenceDiagrams/distances/detail/Munkres.hh>
#include <aleph/persistenceDiagrams/distances/detail/Orthogonal.hh>

#include <algorithm>
#include <iterator>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace aleph
{

namespace detail
{

// Auxiliary structure for describing a pairing between persistence
// diagrams.
struct Pairing
{
  using IndexType = std::size_t;
  using Pair      = std::pair<IndexType, IndexType>;

  bool operator==( const Pairing& other ) const noexcept
  {
    return cost == other.cost;
  }

  double cost;
  std::vector<Pair> pairs;
};

template <
  class DataType,
  class Distance = aleph::distances::InfinityDistance<DataType>
> Pairing optimalPairing( const PersistenceDiagram<DataType>& D1,
                          const PersistenceDiagram<DataType>& D2,
                          DataType power = DataType( 2 ) )
{

  if( D1.dimension() != D2.dimension() )
    throw std::runtime_error( "Dimensions do not coincide" );

  auto size = D1.size() + D2.size();

  distances::detail::Matrix<DataType> costs( size );

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
        d = std::pow( distances::detail::orthogonalDistance<Distance>( p1 ), power );
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
        d = std::pow( distances::detail::orthogonalDistance<Distance>( p1 ), power );
      else
        d = std::numeric_limits<DataType>::max();

      costs( row, col ) = d;

      ++col;
    }

    ++row;
  }

  // Assignment problem solving ----------------------------------------

  distances::detail::Munkres<DataType> solver( costs );

  auto M                                           = solver();
  aleph::math::KahanSummation<DataType> totalCosts = DataType();

  Pairing pairing;

  for( row = IndexType(); row < M.n(); row++ )
  {
    for( col = IndexType(); col < M.n(); col++ )
    {
      if( M( row, col ) == IndexType() )
      {
        // This ensures that pairs are returned in the order dictated by
        // the first persistence diagram.
        pairing.pairs.push_back( std::make_pair( row, col ) );
        totalCosts += costs( row, col );
      }
    }
  }

  pairing.cost = totalCosts;
  return pairing;
}

} // namespace detail

template <class InputIterator> auto mean( InputIterator begin, InputIterator end ) -> typename std::iterator_traits<InputIterator>::value_type
{
  using PersistenceDiagram = typename std::iterator_traits<InputIterator>::value_type;
  using DataType           = typename PersistenceDiagram::DataType;

  std::vector<PersistenceDiagram> persistenceDiagrams( begin, end );
  PersistenceDiagram Y;

  {
    std::random_device rd;
    std::default_random_engine rng( rd() );
    std::uniform_int_distribution<decltype( persistenceDiagrams.size() )> distribution( 0, persistenceDiagrams.size() - 1 );

    Y = persistenceDiagrams.at( distribution( rng ) );
    Y.removeDiagonal();
  }

  aleph::math::KahanSummation<double> cost = 0.0;

  std::vector<detail::Pairing> pairings;
  pairings.resize( decltype(pairings)::size_type( std::distance( begin, end ) ) );

  #pragma omp parallel for
  for( auto it = begin; it < end; ++it )
  {
    auto i      = decltype(pairings)::size_type( std::distance( begin, it ) );
    pairings[i] = detail::optimalPairing( Y, *it );

    #pragma omp critical
    {
      cost += pairings.back().cost;
    }
  }

  bool stop = false;
  while( !stop )
  {
    PersistenceDiagram Z;

    auto k = Y.size();
    for( decltype(k) i = 0; i < k; i++ )
    {
      aleph::math::KahanSummation<DataType> x0 = DataType(); // off-diagonal
      aleph::math::KahanSummation<DataType> y0 = DataType(); // off-diagonal
      aleph::math::KahanSummation<DataType> x1 = DataType();
      aleph::math::KahanSummation<DataType> y1 = DataType();

      using DifferenceType = typename std::iterator_traits<InputIterator>::difference_type;
      DifferenceType index = 0;

      {
        auto point = *std::next( Y.begin(), DifferenceType( i ) );

        // The orthogonal projection is given by 0.5*(x+y). I am using
        // a different calculation to prevent implicit conversions.
        x1 = ( point.x() + point.y() ) / 2;
        y1 = ( point.x() + point.y() ) / 2;
      }

      // Counting the number of off-diagonal points is required in order
      // to ensure that the arithmetical mean is weighted correctly.
      unsigned numOffDiagonalPoints = 0;

      for( auto&& pairing : pairings )
      {
        auto&& diagram = *( std::next( begin, index ) );

        // Off-diagonal assignment
        if( pairing.pairs.at(i).second < diagram.size() )
        {
          auto&& point = *( std::next( diagram.begin(), DifferenceType( pairing.pairs.at(i).second ) ) );

          x0 += point.x();
          y0 += point.y();

          numOffDiagonalPoints++;
        }

        ++index;
      }

      auto x = ( x0 + ( DataType( persistenceDiagrams.size() - numOffDiagonalPoints ) ) * x1 ) / DataType( pairings.size() );
      auto y = ( y0 + ( DataType( persistenceDiagrams.size() - numOffDiagonalPoints ) ) * y1 ) / DataType( pairings.size() );

      Z.add( x,y );
    }

    Y = Z;
    Y.removeDiagonal();

    {
      std::vector<detail::Pairing> newPairings;
      newPairings.resize( decltype(newPairings)::size_type( std::distance( begin, end ) ) );

      aleph::math::KahanSummation<double> newCost = 0.0;

      #pragma omp parallel for
      for( auto it = begin; it < end; ++it )
      {
        auto i         = decltype(newPairings)::size_type( std::distance( begin, it ) );
        newPairings[i] = detail::optimalPairing( Y, *it );

        #pragma omp critical
        {
          newCost += newPairings.back().cost;
        }
      }

      if( newPairings == pairings )
        stop = true;
      else
      {
        pairings.swap( newPairings );
        cost = newCost;
      }
    }
  }

  return Y;
}


} // namespace aleph

#endif
