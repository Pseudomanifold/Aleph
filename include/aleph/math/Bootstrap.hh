#ifndef ALEPH_MATH_BOOTSTRAP_HH__
#define ALEPH_MATH_BOOTSTRAP_HH__

#include <iterator>
#include <random>
#include <vector>

namespace aleph
{

namespace math
{

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
    typedef typename std::iterator_traits<InputIterator>::value_type SampleValueType;
    typedef decltype( functor( begin, end ) ) FunctorValueType;

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
};

} // namespace math

} // namespace aleph

#endif
