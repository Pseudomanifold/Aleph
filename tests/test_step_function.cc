#include "math/StepFunction.hh"

#include "persistenceDiagrams/PersistenceDiagram.hh"
#include "persistenceDiagrams/PersistenceIndicatorFunction.hh"

#include "tests/Base.hh" 

using namespace aleph::math;

template <class T> void testStepFunction()
{
  ALEPH_TEST_BEGIN( "Step function: Basic properties" );

  StepFunction<T> f;
  f.add( 0, 1, 1   );
  f.add( 2, 3, 1   );
  f.add( 3, 4, 2 );

  StepFunction<T> g;
  g.add( 0.5, 0.75, 1 );

  ALEPH_ASSERT_THROW( f(0)   == 1 );
  ALEPH_ASSERT_THROW( f(1)   == 1 );
  ALEPH_ASSERT_THROW( f(1.5) == 0 );
  ALEPH_ASSERT_THROW( f(2)   == 1 );
  ALEPH_ASSERT_THROW( f(3)   == 3 );
  ALEPH_ASSERT_THROW( f(3.5) == 2 );
  ALEPH_ASSERT_THROW( g(0.5) == 1 );
  ALEPH_ASSERT_THROW( g(1.0) == 0 );

  ALEPH_ASSERT_THROW( f.integral() == 4    );
  ALEPH_ASSERT_THROW( g.integral() == 0.25 );

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

  testPersistenceIndicatorFunction<double>();
  testPersistenceIndicatorFunction<float>();
}
