#include <tests/Base.hh>

#include <aleph/math/Bootstrap.hh>

#include <iterator>
#include <numeric>
#include <vector>

auto meanCalculation = [] ( auto begin, auto end )
{
  using T  = typename std::iterator_traits<decltype(begin)>::value_type;
  auto sum = std::accumulate( begin, end, T() );

  return static_cast<double>( sum / std::distance(begin, end) );
};

void testSimple()
{
  std::vector<unsigned> samples = {1,1,1,2,3,4,5,5,5};
  std::vector<double> means;

  auto mean = meanCalculation( samples.begin(), samples.end() );

  ALEPH_ASSERT_EQUAL( mean, 3.0 );

  aleph::math::Bootstrap bootstrap;

  bootstrap.makeReplicates( 200,
                            samples.begin(), samples.end(),
                            meanCalculation,
                            std::back_inserter( means ) );
}

int main(int, char**)
{
  testSimple();
}
