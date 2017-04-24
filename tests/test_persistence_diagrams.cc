#include "tests/Base.hh"

#include "persistenceDiagrams/MultiScaleKernel.hh"
#include "persistenceDiagrams/PersistenceDiagram.hh"

#include <limits>
#include <random>

#include <cmath>

template <class T> aleph::PersistenceDiagram<T> createRandomPersistenceDiagram( unsigned n )
{
  std::random_device rd;
  std::default_random_engine rng( rd() );
  std::uniform_real_distribution<T> distribution( T(0), T( std::nextafter( T(1), std::numeric_limits<T>::max() ) ) );

  using PersistenceDiagram = aleph::PersistenceDiagram<T>;
  PersistenceDiagram D;

  for( unsigned i = 0; i < n; i++ )
  {
    auto x = distribution( rng );
    auto y = distribution( rng );

    if( x > y )
      std::swap( x,y );

    D.add( x,y );
  }

  return D;
}

template <class T> void testMultiScaleKernel()
{
  ALEPH_TEST_BEGIN( "Multi-scale kernel" );

  auto D1 = createRandomPersistenceDiagram<T>( 50 );
  auto D2 = createRandomPersistenceDiagram<T>( 50 );

  auto d1 = aleph::multiScalePseudoMetric( D1, D2, 1.0 );
  auto d2 = aleph::multiScalePseudoMetric( D1, D2, 2.0 );

  ALEPH_ASSERT_THROW( d1 > 0.0 );
  ALEPH_ASSERT_THROW( d2 > 0.0 );

  ALEPH_TEST_END();
}

int main(int, char**)
{
  testMultiScaleKernel<float> ();
  testMultiScaleKernel<double>();
}
