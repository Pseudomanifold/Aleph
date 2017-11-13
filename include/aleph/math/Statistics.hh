#ifndef ALEPH_MATH_STATISTICS_HH__
#define ALEPH_MATH_STATISTICS_HH__

#include <aleph/math/Quantiles.hh>

#include <iterator>
#include <limits>

#include <cmath>

namespace aleph
{

namespace math
{

/**
  Given a range of numbers, calculates the sample variance of the data.

  @param begin Input iterator to begin of range
  @param end   Input iterator to end of range

  @returns Sample variance
*/

template <class InputIterator> double sampleVariance( InputIterator begin,
                                                      InputIterator end)
{
  if(    begin == end
      || std::distance( begin, end ) == 1  )
  {
    return( std::numeric_limits<double>::quiet_NaN() );
  }

  double mean = math::mean( begin, end );
  double sum  = 0.0;

  for( InputIterator it = begin; it != end; ++it )
  {
    double d  = static_cast<double>( *it );
    sum      += ( d - mean ) * ( d - mean );
  }

  return ( 1.0 / static_cast<double>( std::distance( begin, end ) - 1 ) ) * sum;
}

/**
  Given a range of numbers, calculates the sample standard deviation of the
  data.

  @param begin Input iterator to begin of range
  @param end   Input iterator to end of range

  @returns Sample standard deviation
*/

template <class InputIterator> double sampleStandardDeviation( InputIterator begin,
                                                               InputIterator end )
{
  return std::sqrt( sampleVariance( begin, end ) );
}

} // namespace math

} // namespace aleph

#endif
