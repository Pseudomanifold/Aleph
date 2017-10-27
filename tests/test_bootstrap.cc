#include <tests/Base.hh>

#include <aleph/math/Bootstrap.hh>

#include <iterator>
#include <numeric>
#include <vector>

auto meanCalculation = [] ( auto begin, auto end )
{
  using T  = typename std::iterator_traits<decltype(begin)>::value_type;
  auto sum = std::accumulate( begin, end, T() );

  return static_cast<double>( sum / static_cast<double>( std::distance(begin, end) ) );
};

void testSimple()
{
  ALEPH_TEST_BEGIN( "Bootstrap: Confidence intervals" );

  // Data and some of the estimates were taken from an MIT course [1],
  // even though their methodology is slightly different.
  //
  // [1]: https://ocw.mit.edu/courses/mathematics/18-05-introduction-to-probability-and-statistics-spring-2014/readings/MIT18_05S14_Reading24.pdf

  std::vector<unsigned> samples = {30,37,36,43,42,43,43,46,41,42};
  std::vector<double> means;

  auto mean                    = meanCalculation( samples.begin(), samples.end() );
  unsigned numBootstrapSamples = 1000;

  ALEPH_ASSERT_EQUAL( mean, 40.3 );

  aleph::math::Bootstrap bootstrap;

  bootstrap.makeReplicates( numBootstrapSamples,
                            samples.begin(), samples.end(),
                            meanCalculation,
                            std::back_inserter( means ) );

  ALEPH_ASSERT_EQUAL( means.size(), numBootstrapSamples );

  // Checking the basic confidence interval of the sample --------------
  //
  // This indicates that the basic confidence interval is not given very
  // specific information (at least not for these data).

  auto basicConfidenceInterval
    = bootstrap.basicConfidenceInterval( numBootstrapSamples,
                                         0.20,
                                         samples.begin(), samples.end(),
                                         meanCalculation );

  ALEPH_ASSERT_THROW( basicConfidenceInterval.first  >= 37.0 );
  ALEPH_ASSERT_THROW( basicConfidenceInterval.second <= 43.0 );

  // Checking the percentile confidence interval -----------------------
  //
  // For laziness reasons, I am using the same upper and lower bounds as
  // for the basic confidence interval.
  //
  // The percentile confidence interval works better when the sample
  // size is larger.

  auto percentileConfidenceInterval
    = bootstrap.percentileConfidenceInterval( numBootstrapSamples,
                                              0.20,
                                              samples.begin(), samples.end(),
                                              meanCalculation );

  ALEPH_ASSERT_THROW( percentileConfidenceInterval.first  >= 37.0 );
  ALEPH_ASSERT_THROW( percentileConfidenceInterval.second <= 43.0 );

  // Checking the student-t confidence interval ------------------------

  auto studentConfidenceInterval
    = bootstrap.studentConfidenceInterval( numBootstrapSamples,
                                           0.20,
                                           samples.begin(), samples.end(),
                                           meanCalculation );

  ALEPH_ASSERT_THROW( studentConfidenceInterval.first  >= 37.0 );
  ALEPH_ASSERT_THROW( studentConfidenceInterval.second <= 43.0 );

  ALEPH_TEST_END();
}

void testStandardError()
{
  ALEPH_TEST_BEGIN( "Bootstrap: Standard error" );

  std::vector<unsigned short> samples = {
    61, 88, 89, 89, 90, 92, 93, 94, 98, 98, 101, 102, 105, 108,
    109, 113, 114, 115, 120, 138
  };

  std::vector<double> means;

  unsigned numBootstrapSamples = 10000;

  aleph::math::Bootstrap bootstrap;

  auto se = bootstrap.standardError(
    numBootstrapSamples,
    samples.begin(), samples.end(),
    meanCalculation );

  ALEPH_ASSERT_THROW( std::abs( se - 3.50 ) < 0.5 );

  ALEPH_TEST_END();
}

int main(int, char**)
{
  testSimple();
  testStandardError();
}
