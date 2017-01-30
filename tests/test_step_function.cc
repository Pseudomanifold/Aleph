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

  ALEPH_ASSERT_THROW( f(0)   == 1 );
  ALEPH_ASSERT_THROW( f(1)   == 1 );
  ALEPH_ASSERT_THROW( f(1.5) == 0 );
  ALEPH_ASSERT_THROW( f(2)   == 1 );
  ALEPH_ASSERT_THROW( f(3)   == 2 );
  ALEPH_ASSERT_THROW( f(3.5) == 2 );
  ALEPH_ASSERT_THROW( g(0.5) == 1 );
  ALEPH_ASSERT_THROW( g(1.0) == 0 );

  ALEPH_ASSERT_THROW( f.integral() == T(4.00) );
  ALEPH_ASSERT_THROW( g.integral() == T(0.25) );

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

  StepFunction<T> f;
  f.add( 0, 1, 1 );

  StepFunction<T> g;
  g.add( 0.25, 0.75, 1 );

  auto h = f+g;

  ALEPH_ASSERT_THROW( h( T(0   ) ) == 1 );
  ALEPH_ASSERT_THROW( h( T(0.20) ) == 1 );
  ALEPH_ASSERT_THROW( h( T(0.25) ) == 2 );
  ALEPH_ASSERT_THROW( h( T(0.50) ) == 2 );
  ALEPH_ASSERT_THROW( h( T(0.75) ) == 2 );
  ALEPH_ASSERT_THROW( h( T(0.80) ) == 1 );
  ALEPH_ASSERT_THROW( h( T(1.00) ) == 1 );

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

  auto f = aleph::persistenceIndicatorFunction( pd1);
  auto g = aleph::persistenceIndicatorFunction( pd2 );

  ALEPH_TEST_END();
}

int main()
{
  testStepFunction<double>();
  testStepFunction<float> ();

  testStepFunctionAddition<double>();
  testStepFunctionAddition<float> ();

  testPersistenceIndicatorFunction<double>();
  testPersistenceIndicatorFunction<float>();
}
