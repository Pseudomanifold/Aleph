#ifndef ALEPH_MATH_BOOTSTRAP_HH__
#define ALEPH_MATH_BOOTSTRAP_HH__

#include <aleph/math/KahanSummation.hh>

#include <boost/math/distributions/students_t.hpp>

#include <algorithm>
#include <iterator>
#include <random>
#include <vector>

#include <cmath>

namespace aleph
{

namespace math
{

/**
  @class Bootstrap
  @brief Generic bootstrap functor

  This functor provides a generic interface for performing bootstrap
  operations on *arbitrary* data, using an *arbitrary* statistic for
  testing. Several convenience functions for estimating *confidence*
  values are provided.
*/

class Bootstrap
{
public:

  /**
    Given a range of data of some type, calculates a set of bootstrap replicates
    for a desired statistic. This function will not perform any type conversions
    in order to preserve all original types. The type of the output data depends
    on the return value type of the functor.

    @param[in]  numSamples samples Number of bootstrap samples
    @param[in]  begin      Input iterator to begin of data range
    @param[in]  end        Input iterator to end of data range
    @param[in]  functor    Functor for calculating a statistic on the replicate
    @param[out] result     Output iterator for storing the results
  */

  template <class InputIterator, class OutputIterator, class Functor>
  void makeReplicates( unsigned numSamples,
                       InputIterator begin, InputIterator end,
                       Functor functor,
                       OutputIterator result )
  {
    using SampleValueType  = typename std::iterator_traits<InputIterator>::value_type;
    using FunctorValueType = decltype( functor( begin, end ) );

    std::vector<SampleValueType> samples( begin, end );

    std::random_device rd;
    std::mt19937 rng( rd() );

    std::uniform_int_distribution<std::size_t> distribution( 0, samples.size() - 1 );

    std::vector<FunctorValueType> replicates;
    replicates.reserve( numSamples );

    for( unsigned sampleIndex = 0; sampleIndex < numSamples; sampleIndex++ )
    {
      std::vector<SampleValueType> sample;
      sample.reserve( samples.size() );

      for( std::size_t i = 0; i < samples.size(); i++ )
        sample.push_back( samples.at( distribution( rng ) ) );

      replicates.push_back( functor( sample.begin(),
                                     sample.end() ) );
    }

    std::copy( replicates.begin(), replicates.end(), result );
  }

  /**
    Calculates a bootstrap estimate of the standard error of a test
    statistics on a data set.

    @param numSamples Number of bootstrap samples
    @param begin      Input iterator to begin of data range
    @param end        Input iterator to end of data range
    @param functor    Functor to describe the test statistics that is to be
                      calculated on the data range. The functor has to take
                      the `value_type` of the range as a parameter, because
                      the function will show the individual samples to it.\n

                      A good example of a functor is the following *mean*
                      functor:\n

    \code{.cpp}
    auto meanCalculation = [] ( auto begin, auto end )
    {
      using T  = typename std::iterator_traits<decltype(begin)>::value_type;
      auto sum = std::accumulate( begin, end, T() );

      return static_cast<double>( sum / static_cast<double>( std::distance(begin, end) ) );
    };
    \endcode
                     Note that the functor requires `C++14` because the use
                     of `auto` in a lambda expression.

    @returns Bootstrapped estimated of the standard error
  */

  template <class InputIterator, class Functor>
  double standardError( unsigned numSamples,
                        InputIterator begin, InputIterator end,
                        Functor functor )
  {
    using FunctorValueType = decltype( functor(begin, end) );

    std::vector<FunctorValueType> estimates;
    estimates.reserve( numSamples );

    this->makeReplicates( numSamples,
                          begin, end,
                          functor,
                          std::back_inserter( estimates ) );

    using namespace aleph::math;

    double mean  = static_cast<double>( accumulate_kahan( estimates.begin(), estimates.end(), FunctorValueType() ) );
    mean        /= numSamples;

    std::vector<double> sampleSquaredDeviations;
    sampleSquaredDeviations.reserve( numSamples );

    for( auto&& estimate : estimates )
    {
      auto delta  = mean - estimate;
      delta      *= delta;

      sampleSquaredDeviations.push_back( delta );
    }

    double sumSquaredDeviations
      = accumulate_kahan_sorted(
          sampleSquaredDeviations.begin(),
          sampleSquaredDeviations.end(),
          0.0 );

    sumSquaredDeviations /= (numSamples - 1);
    return std::sqrt( static_cast<double>( sumSquaredDeviations ) );
  }

  template <class InputIterator, class Functor>
  auto basicConfidenceInterval( unsigned numSamples,
                                double alpha,
                                InputIterator begin, InputIterator end,
                                Functor functor ) -> std::pair< decltype( functor(begin, end) ), decltype( functor(begin, end) ) >
  {
    auto theta             = functor( begin, end );
    using FunctorValueType = decltype( theta );

    std::vector<FunctorValueType> estimates;
    estimates.reserve( numSamples );

    this->makeReplicates( numSamples,
                          begin, end,
                          functor,
                          std::back_inserter( estimates ) );

    std::sort( estimates.begin(), estimates.end() );

    auto upperPercentile = alpha / 2;
    auto upperEstimate   = estimates.at( Bootstrap::index( numSamples, upperPercentile ) );
    auto lowerPercentile = 1 - upperPercentile;
    auto lowerEstimate   = estimates.at( Bootstrap::index( numSamples, lowerPercentile ) );

    return std::make_pair( 2*theta - lowerEstimate, 2*theta - upperEstimate );
  }

  template <class InputIterator, class Functor>
  auto percentileConfidenceInterval( unsigned numSamples,
                                     double alpha,
                                     InputIterator begin, InputIterator end,
                                     Functor functor ) -> std::pair< decltype( functor(begin, end) ), decltype( functor(begin, end) ) >
  {
    using FunctorValueType = decltype( functor( begin, end ) );

    std::vector<FunctorValueType> estimates;
    estimates.reserve( numSamples );

    this->makeReplicates( numSamples,
                          begin, end,
                          functor,
                          std::back_inserter( estimates ) );

    std::sort( estimates.begin(), estimates.end() );

    auto lowerPercentile = alpha / 2;
    auto lowerEstimate   = estimates.at( Bootstrap::index( numSamples, lowerPercentile ) );
    auto upperPercentile = 1 - lowerPercentile;
    auto upperEstimate   = estimates.at( Bootstrap::index( numSamples, upperPercentile ) );

    return std::make_pair( lowerEstimate, upperEstimate );
  }

  template <class InputIterator, class Functor>
  auto studentConfidenceInterval( unsigned numSamples,
                                  double alpha,
                                  InputIterator begin, InputIterator end,
                                  Functor functor ) -> std::pair< decltype( functor(begin, end) ), decltype( functor(begin, end) ) >
  {
    auto theta = functor( begin, end );

    boost::math::students_t distribution( static_cast<double>( std::distance( begin, end ) - 1 ) );

    auto tl = boost::math::quantile( distribution, 1 - alpha );
    auto tu = boost::math::quantile( distribution,     alpha );
    auto se = this->standardError( numSamples,
                                   begin, end,
                                   functor );

    return std::make_pair( theta - tl * se, theta - tu * se );
  }

private:

  /** Calculates index at a certain percentile of the data */
  static unsigned index( unsigned int samples, double alpha )
  {
    // This accounts for rounding and works regardless of whether
    // the product samples * alpha is an integer or not. Note the
    // offset of -1. It is required because, say, the 100th value
    // is at index 99 of the vector.
    return static_cast<unsigned>( samples * alpha + 0.5 ) - 1;
  }
};

} // namespace math

} // namespace aleph

#endif
