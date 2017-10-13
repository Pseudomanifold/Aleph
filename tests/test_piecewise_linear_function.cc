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

  ALEPH_ASSERT_EQUAL( f(0), T(0) );
  ALEPH_ASSERT_EQUAL( f(1), T(1) );
  ALEPH_ASSERT_EQUAL( g(0), T(0) );

  ALEPH_ASSERT_EQUAL( f( T(0.5) ),  0.5 );
  ALEPH_ASSERT_EQUAL( g( T(0.5) ), -1   );

  ALEPH_ASSERT_THROW( g( T(0.2) ) <  0 );
  ALEPH_ASSERT_THROW( g( T(0.2) ) > -1 );

  ALEPH_TEST_END();
}

int main()
{
  testBasicProperties<float> ();
  testBasicProperties<double>();
}
