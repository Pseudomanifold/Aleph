#include <tests/Base.hh>

#include <aleph/math/PiecewiseLinearFunction.hh>

#include <iterator>
#include <set>
#include <vector>

#include <cmath>

using namespace aleph::math;

template <class T> void testBasic()
{
  ALEPH_TEST_BEGIN( "Piecewise linear function: Basic properties & operations" );

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

  auto h = f+g;

  ALEPH_ASSERT_EQUAL( h(0),   T(0) );
  ALEPH_ASSERT_EQUAL( h(1),   T(1) );
  ALEPH_ASSERT_EQUAL( h(0.5), T(-0.5) );

  h -= g;

  ALEPH_ASSERT_EQUAL( h(0), T(0) );
  ALEPH_ASSERT_EQUAL( h(1), T(1) );
  ALEPH_ASSERT_EQUAL( h( T(0.5) ),  0.5 );

  h *= 2;

  ALEPH_ASSERT_EQUAL( h(1), 2*T(1) );
  ALEPH_ASSERT_EQUAL( h( T(0.5) ), 2*0.5 );

  h /= 2;

  ALEPH_ASSERT_EQUAL( h(0), T(0) );
  ALEPH_ASSERT_EQUAL( h(1), T(1) );
  ALEPH_ASSERT_EQUAL( h( T(0.5) ),  0.5 );

  h = -h;

  ALEPH_ASSERT_EQUAL( h(0), T(0) );
  ALEPH_ASSERT_EQUAL( h(1), -T(1) );
  ALEPH_ASSERT_EQUAL( h( T(0.5) ), -0.5 );

  ALEPH_ASSERT_THROW( h != -h );
  ALEPH_ASSERT_THROW( h == -f );

  ALEPH_ASSERT_THROW( h.abs() == f );

  ALEPH_ASSERT_EQUAL( f.integral(), T(0.5) );
  ALEPH_ASSERT_EQUAL( f.integral(), (-f).integral() );
  ALEPH_ASSERT_EQUAL( g.integral(), T(0.25) );
  ALEPH_ASSERT_EQUAL( f.integral(), h.integral() );

  ALEPH_ASSERT_THROW( (f+g).integral() < (f.integral() + g.integral() ) );

  {
    auto f = PiecewiseLinearFunction<T>();

    ALEPH_ASSERT_THROW( f == PiecewiseLinearFunction<T>() );
    ALEPH_ASSERT_THROW( f+f == f );
    ALEPH_ASSERT_EQUAL( f.integral(), PiecewiseLinearFunction<T>().integral() );
  }

  ALEPH_TEST_END();
}

int main()
{
  testBasic<float> ();
  testBasic<double>();
}
