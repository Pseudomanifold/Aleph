#include "math/StepFunction.hh"

#include "persistenceDiagrams/PersistenceDiagram.hh"
#include "persistenceDiagrams/PersistenceIndicatorFunction.hh"

#include "tests/Base.hh"

#include <set>

#include <cmath>

using namespace aleph::math;

template <class T> bool almostEqual( T x, T y )
{
  auto difference = std::abs( x-y );
  x               = std::abs( x );
  y               = std::abs( y );

  auto largest    = x < y ? y : x;

  return difference <= largest * 2 * std::numeric_limits<T>::epsilon();
}

template <class T> void testStepFunction()
{
  ALEPH_TEST_BEGIN( "Step function: Basic properties" );

  StepFunction<T> f;
  f.add( 0, 1, 1 );
  f.add( 2, 3, 1 );
  f.add( 3, 4, 2 );

  StepFunction<T> g;
  g.add( 0.5, 0.75, 1 );

  ALEPH_ASSERT_EQUAL( f(0)  , 1 );
  ALEPH_ASSERT_EQUAL( f(1)  , 1 );
  ALEPH_ASSERT_EQUAL( f(1.5), 0 );
  ALEPH_ASSERT_EQUAL( f(2)  , 1 );
  ALEPH_ASSERT_EQUAL( f(3)  , 2 );
  ALEPH_ASSERT_EQUAL( f(3.5), 2 );
  ALEPH_ASSERT_EQUAL( f(4.0), 2);

  ALEPH_ASSERT_EQUAL( g(0.5),1 );
  ALEPH_ASSERT_EQUAL( g(1.0), 0 );

  ALEPH_ASSERT_EQUAL( f.integral(), T(4.00) );
  ALEPH_ASSERT_EQUAL( g.integral(), T(0.25) );

  auto h = f+g;

  std::set<T> D;
  f.domain( std::inserter( D, D.begin() ) );
  g.domain( std::inserter( D, D.begin() ) );

  for( auto&& x : D )
    ALEPH_ASSERT_THROW( h(x) != 0 );

  ALEPH_TEST_END();
}

template <class T> void testStepFunctionAddition()
{
  ALEPH_TEST_BEGIN( "Step function: Addition" );

  // Case 1: Inclusion -------------------------------------------------
  //
  // Use indicator intervals that are included in each other without any
  // overlaps.
  {
    StepFunction<T> f;
    f.add( 0, 1, 1 );

    StepFunction<T> g;
    g.add( 0.25, 0.75, 2 );

    auto h = f+g;

    ALEPH_ASSERT_EQUAL( h( T(0   ) ), 1 );
    ALEPH_ASSERT_EQUAL( h( T(0.20) ), 1 );
    ALEPH_ASSERT_EQUAL( h( T(0.25) ), 3 );
    ALEPH_ASSERT_EQUAL( h( T(0.50) ), 3 );
    ALEPH_ASSERT_EQUAL( h( T(0.75) ), 3 );
    ALEPH_ASSERT_EQUAL( h( T(0.80) ), 1 );
    ALEPH_ASSERT_EQUAL( h( T(1.00) ), 1 );
  }

  // Case 2: Overlap ---------------------------------------------------
  //
  // Use indicator intervals that overlap without being equal. This must
  // only result in updates within the 'shared' region of the functions.
  {
    StepFunction<T> f;
    f.add( 0, 1, 1 );

    StepFunction<T> g;
    g.add( 0.50, 1.50, 2 );

    auto h = f+g;

    ALEPH_ASSERT_EQUAL( h( T(0   ) ), 1 );
    ALEPH_ASSERT_EQUAL( h( T(0.40) ), 1 );
    ALEPH_ASSERT_EQUAL( h( T(0.50) ), 3 );
    ALEPH_ASSERT_EQUAL( h( T(0.75) ), 3 );
    ALEPH_ASSERT_EQUAL( h( T(1.00) ), 3 );
    ALEPH_ASSERT_EQUAL( h( T(1.10) ), 2 );
    ALEPH_ASSERT_EQUAL( h( T(1.50) ), 2 );
  }

  // Case 3: Touching ---------------------------------------------------
  //
  // Use indicator intervals whose intervals touch. This is interesting
  // insofar it requires creating a new interval directly subsequent to
  // the critical point
  {
    StepFunction<T> f;
    f.add( 0, 1, 1 );

    StepFunction<T> g;
    g.add( 1, 2, 2 );

    auto h = f+g;

    ALEPH_ASSERT_EQUAL( h( T(0   ) ), 1 );
    ALEPH_ASSERT_EQUAL( h( T(0.50) ), 1 );
    ALEPH_ASSERT_EQUAL( h( T(1.00) ), 3 );
    ALEPH_ASSERT_EQUAL( h( T(1.01) ), 2 );
    ALEPH_ASSERT_EQUAL( h( T(1.50) ), 2 );
    ALEPH_ASSERT_EQUAL( h( T(2.00) ), 2 );
  }

  // Case 4: Equality --------------------------------------------------
  //
  // If the functions fully coincide, this should be equivalent to scalar
  // multiplication.

  {
    StepFunction<T> f;
    f.add( 1,2, 1);

    auto g = f*2;
    auto h = f+f;

    ALEPH_ASSERT_EQUAL( h(0  ), g(0  ) );
    ALEPH_ASSERT_EQUAL( h(1  ), g(1  ) );
    ALEPH_ASSERT_EQUAL( h(1.5), g(1.5) );
    ALEPH_ASSERT_EQUAL( h(2  ), g(2  ) );
  }

  ALEPH_TEST_END();
}

template <class T> void testStepFunctionNegation()
{
  ALEPH_TEST_BEGIN( "Step function: Negation" );

  StepFunction<T> f;
  f.add( 0, 1, 1 );
  f.add( 2, 3, 1 );
  f.add( 3, 4, 2 );

  auto g =  f * -1;
  auto h = -f;

  ALEPH_ASSERT_EQUAL( g(0)  , -1 );
  ALEPH_ASSERT_EQUAL( g(1)  , -1 );
  ALEPH_ASSERT_EQUAL( g(1.5), -0 );
  ALEPH_ASSERT_EQUAL( g(2)  , -1 );
  ALEPH_ASSERT_EQUAL( g(3)  , -2 );
  ALEPH_ASSERT_EQUAL( g(3.5), -2 );
  ALEPH_ASSERT_EQUAL( g(4.0), -2);

  ALEPH_ASSERT_EQUAL( g(0)  , h(0)   );
  ALEPH_ASSERT_EQUAL( g(1)  , h(1)   );
  ALEPH_ASSERT_EQUAL( g(1.5), h(1.5) );
  ALEPH_ASSERT_EQUAL( g(2)  , h(2)   );
  ALEPH_ASSERT_EQUAL( g(3)  , h(3)   );
  ALEPH_ASSERT_EQUAL( g(3.5), h(3.5) );
  ALEPH_ASSERT_EQUAL( g(4.0), h(4.0) );

  ALEPH_TEST_END();
}

template <class T> void testStepFunctionNormalization()
{
  ALEPH_TEST_BEGIN( "Step function: Normalization" );

  StepFunction<T> f;
  f.add( 0, 1, 1 );
  f.add( 2, 3, 1 );
  f.add( 3, 4, 2 );

  auto g = normalize( f );

  ALEPH_ASSERT_EQUAL( f(0)  , 1 );
  ALEPH_ASSERT_EQUAL( f(1)  , 1 );
  ALEPH_ASSERT_EQUAL( f(1.5), 0 );
  ALEPH_ASSERT_EQUAL( f(2)  , 1 );
  ALEPH_ASSERT_EQUAL( f(3)  , 2 );
  ALEPH_ASSERT_EQUAL( f(3.5), 2 );
  ALEPH_ASSERT_EQUAL( f(4.0), 2);

  ALEPH_ASSERT_EQUAL( g(0)  , 0.5 );
  ALEPH_ASSERT_EQUAL( g(1)  , 0.5 );
  ALEPH_ASSERT_EQUAL( g(1.5), 0.0 );
  ALEPH_ASSERT_EQUAL( g(2)  , 0.5 );
  ALEPH_ASSERT_EQUAL( g(3)  , 1.0 );
  ALEPH_ASSERT_EQUAL( g(3.5), 1.0 );
  ALEPH_ASSERT_EQUAL( g(4.0), 1.0);

  ALEPH_TEST_END();
}

template <class T> void testPersistenceIndicatorFunction()
{
  using PersistenceDiagram = aleph::PersistenceDiagram<T>;

  ALEPH_TEST_BEGIN( "Persistence indicator function" );

  PersistenceDiagram pd1;
  pd1.add(1,   2  );
  pd1.add(1.5, 2.5);
  pd1.add(2,   3  );

  PersistenceDiagram pd2;
  pd2.add(1,   2  );
  pd2.add(3,   4  );

  PersistenceDiagram pd3;
  pd3.add(0,1);
  pd3.add(1,2);

  PersistenceDiagram pd4;
  pd4.add(0,1);
  pd4.add(0,6);
  pd4.add(1,2);
  pd4.add(2,3);
  pd4.add(3,6);
  pd4.add(5,8);

  auto f = aleph::persistenceIndicatorFunction( pd1);
  auto g = aleph::persistenceIndicatorFunction( pd2 );
  auto h = aleph::persistenceIndicatorFunction( pd3 );
  auto i = aleph::persistenceIndicatorFunction( pd4 );

  ALEPH_ASSERT_EQUAL( h(0), 1 );
  ALEPH_ASSERT_EQUAL( h(2), 1 );

  ALEPH_ASSERT_EQUAL( i(T(0.1)), 2 );
  ALEPH_ASSERT_EQUAL( i(T(1.1)), 2 );
  ALEPH_ASSERT_EQUAL( i(T(2.1)), 2 );
  ALEPH_ASSERT_EQUAL( i(T(3.1)), 2 );
  ALEPH_ASSERT_EQUAL( i(T(5.1)), 3 );
  ALEPH_ASSERT_EQUAL( i(T(6.1)), 1 );
  ALEPH_ASSERT_EQUAL( i(T(8.1)), 0 );

  ALEPH_TEST_END();
}

int main()
{
  testStepFunction<double>();
  testStepFunction<float> ();

  testStepFunctionAddition<double>();
  testStepFunctionAddition<float> ();

  testStepFunctionNegation<double>();
  testStepFunctionNegation<float> ();

  testStepFunctionNormalization<double>();
  testStepFunctionNormalization<float> ();

  testPersistenceIndicatorFunction<double>();
  testPersistenceIndicatorFunction<float>();
}
