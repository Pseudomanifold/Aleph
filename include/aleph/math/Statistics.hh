#ifndef ALEPH_MATH_STATISTICS_HH__
#define ALEPH_MATH_STATISTICS_HH__

#include <iterator>
#include <limits>
#include <numeric>

#include <cmath>

namespace aleph
{

namespace math
{

/**
  Given a range of numbers, calculates the sample mean of the data.

  @param begin Input iterator to begin of range
  @param end   Input iterator to end of range

  @returns Sample mean
*/

template <class InputIterator> double sampleMean( InputIterator begin,
                                                  InputIterator end )
{
  if( begin == end )
    return std::numeric_limits<double>::quiet_NaN();

  return std::accumulate( begin, end, 0.0 ) / static_cast<double>( std::distance( begin, end ) );
}

/**
  Given two ranges of numbers, calculates their sample covariance and
  returns it.

  @param begin Input iterator to begin of range
  @param end   Input iterator to end of range

  @returns Sample covariance
*/

template <
  class InputIterator1,
  class InputIterator2
> double sampleCovariance( InputIterator1 begin1, InputIterator1 end1,
                           InputIterator2 begin2, InputIterator2 end2 )
{
  if(    begin1 == end1
      || begin2 == end2
      || std::distance( begin1, end1 ) == 1
      || std::distance( begin2, end2 ) == 1  )
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double mean1 = sampleMean( begin1, end1 );
  double mean2 = sampleMean( begin2, end2 );
  double sum   = 0.0;

  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;

  for( ; it1 != end1 && it2 != end2; ++it1, ++it2 )
  {
    double d1  = static_cast<double>( *it1 ) - mean1;
    double d2  = static_cast<double>( *it2 ) - mean2;
    sum       += d1 * d2;
  }

  return sum;
}

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
    return std::numeric_limits<double>::quiet_NaN();
  }

  double mean = sampleMean( begin, end );
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

/**
  Given a range of numbers, calculates the sample skewness of the data.

  @param begin Input iterator to begin of range
  @param end   Input iterator to end of range

  @returns Sample skewness
*/

template <class InputIterator> double sampleSkewness( InputIterator begin,
                                                      InputIterator end )
{
  if( begin == end )
    std::numeric_limits<double>::quiet_NaN();

  double mean        = sampleMean( begin, end );
  double numerator   = 0.0;
  double denominator = 0.0;

  for( InputIterator it = begin; it != end; ++it )
  {
    numerator   += std::pow( static_cast<double>( *it - mean ), 3.0 );
    denominator += ( *it - mean ) * ( *it - mean );
  }

  numerator   /= static_cast<double>( std::distance( begin, end ) );

  denominator /= static_cast<double>( std::distance( begin, end ) );
  denominator  = std::pow( denominator, 3.0 / 2.0 );

  return numerator / denominator;
}

/**
  Given a range of numbers, calculates the sample excess kurtosis of the
  data.

  @param begin Input iterator to begin of data
  @param end   Input iterator to end of data

  @returns Sample excess kurtosis
*/

template <class InputIterator> double sampleKurtosis( InputIterator begin,
                                                      InputIterator end )
{
  if( begin == end )
    return std::numeric_limits<double>::quiet_NaN();

  double mean = sampleMean( begin, end );
  double m4   = 0.0;
  double m2   = 0.0;

  double n    = static_cast<double>( std::distance( begin, end ) );

  for( InputIterator it = begin; it != end; ++it )
  {
    m4 += std::pow( static_cast<double>( *it - mean ), 4.0 );
    m2 += std::pow( static_cast<double>( *it - mean ), 2.0 );
  }

  m4 /= n;
  m2 /= n;

  return m4 / (m2 * m2 ) - 3.0;
}

} // namespace math

} // namespace aleph

#endif
