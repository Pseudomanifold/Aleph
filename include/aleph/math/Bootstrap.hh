#ifndef ALEPH_MATH_BOOTSTRAP_HH__
#define ALEPH_MATH_BOOTSTRAP_HH__

#include <aleph/math/KahanSummation.hh>

#include <algorithm>
#include <iterator>
#include <random>
#include <vector>

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
