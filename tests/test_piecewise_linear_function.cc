#include <tests/Base.hh>

#include <aleph/math/PiecewiseLinearFunction.hh>

#include <iterator>
#include <set>
#include <vector>

#include <cmath>

using namespace aleph::math;

template <class T> void testBasicProperties()
{
  ALEPH_TEST_BEGIN( "Piecewise linear function: Basic properties" );

  std::vector< std::pair<T, T> > points1 = {
    {0,0},
    {1,1}
  };

  std::vector< std::pair<T, T> > points2 = {
    {0,   0},
    {0.5,-1}
  };

  PiecewiseLinearFunction<T> f( points1.begin(), points1.end() );
  PiecewiseLinearFunction<T> g( points2.begin(), points2.end() );

  ALEPH_TEST_END();
}

int main()
{
  testBasicProperties<float> ();
  testBasicProperties<double>();
}
