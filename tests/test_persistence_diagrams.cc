#include <tests/Base.hh>

#include <aleph/persistenceDiagrams/Mean.hh>
#include <aleph/persistenceDiagrams/Norms.hh>
#include <aleph/persistenceDiagrams/MultiScaleKernel.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistenceDiagrams/distances/NearestNeighbour.hh>
#include <aleph/persistenceDiagrams/distances/Wasserstein.hh>

#include <limits>
#include <random>
#include <vector>

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

template <class T> void testFrechetMean()
{
  using PersistenceDiagram = aleph::PersistenceDiagram<T>;

  ALEPH_TEST_BEGIN( "Persistence diagram mean");

  unsigned n = 10;

  std::vector<PersistenceDiagram> diagrams;
  diagrams.reserve( n );

  for( decltype(n) i = 0; i < n; i++ )
    diagrams.emplace_back( createRandomPersistenceDiagram<T>( 25 ) );

  auto D = aleph::mean( diagrams.begin(), diagrams.end() );
  auto P = aleph::totalPersistence( D );
  auto p = std::sqrt(25.0) * std::sqrt(0.50); // simple estimate of the mean value
                                              // for the total persistence

  ALEPH_ASSERT_THROW( D.size() > 0 );
  ALEPH_ASSERT_THROW( std::abs( P - p ) < 2.0 );

  ALEPH_TEST_END();
}

template <class T> void testMultiScaleKernel()
{
  ALEPH_TEST_BEGIN( "Multi-scale kernel" );

  auto D1 = createRandomPersistenceDiagram<T>( 50 );
  auto D2 = createRandomPersistenceDiagram<T>( 50 );

  auto d1 = aleph::multiScalePseudoMetric( D1, D2, 1.0 );
  auto d2 = aleph::multiScalePseudoMetric( D1, D2, 2.0 );
  auto d3 = aleph::distances::wassersteinDistance( D1, D2, T(1) );

  // Non-negativity
  ALEPH_ASSERT_THROW( d1 > 0.0 );
  ALEPH_ASSERT_THROW( d2 > 0.0 );
  ALEPH_ASSERT_THROW( d3 > 0.0 );

  // Symmetry
  ALEPH_ASSERT_EQUAL( aleph::multiScalePseudoMetric(D1, D1, 1.0), aleph::multiScalePseudoMetric(D2, D2, 1.0) );
  ALEPH_ASSERT_EQUAL( aleph::multiScalePseudoMetric(D1, D1, 1.0), 0.0 );

  // Stability
  ALEPH_ASSERT_THROW( d1 < 1.0 / ( 1.0 / ( 1.0 * std::sqrt( 8.0 * M_PI ) ) * d3 ) );
  ALEPH_ASSERT_THROW( d2 < 1.0 / ( 1.0 / ( 2.0 * std::sqrt( 8.0 * M_PI ) ) * d3 ) );

  ALEPH_TEST_END();
}

template <class T> void testNearestNeighbourDistance()
{
  ALEPH_TEST_BEGIN( "Nearest neighbour distance" );

  auto D1 = createRandomPersistenceDiagram<T>( 50 );
  auto D2 = createRandomPersistenceDiagram<T>( 50 );

  auto dNearestNeighbour = aleph::distances::nearestNeighbourDistance( D1, D2 );
  auto dWasserstein      = aleph::distances::wassersteinDistance( D1, D2, T(1) );

  ALEPH_ASSERT_THROW( dNearestNeighbour < dWasserstein );

  ALEPH_TEST_END();
}


template <class T> void testWassersteinDistance()
{
  using Diagram = aleph::PersistenceDiagram<T>;
  using namespace aleph::distances;

  Diagram D1;
  D1.add( T(0.9), T(1.0) );
  D1.add( T(1.9), T(2.0) );
  D1.add( T(2.9), T(3.0) );
  D1.add( T(3.9), T(4.0) );

  {
    auto d11 = wassersteinDistance( D1, D1 );

    assert( d11 >= T() );
    assert( d11 == T() );

    ALEPH_ASSERT_THROW( d11 >= T() );
    ALEPH_ASSERT_THROW( d11 == T() );
  }

  Diagram D2;
  D2.add( T(0.9), T(1.0) );
  D2.add( T(1.9), T(2.0) );
  D2.add( T(2.9), T(3.0) );
  D2.add( T(3.9), T(9.9) );

  {
    auto d12 = wassersteinDistance( D1, D2 );
    auto d21 = wassersteinDistance( D2, D1 );

    ALEPH_ASSERT_THROW( d12 > T() );
    ALEPH_ASSERT_THROW( d21 > T() );

    ALEPH_ASSERT_THROW( d12 == d21 );
    ALEPH_ASSERT_THROW( std::abs( d12 -  T( 3.05 ) ) < 1e-8 );
  }
}

int main(int, char**)
{
  testFrechetMean<float> ();
  testFrechetMean<double>();

  testMultiScaleKernel<float> ();
  testMultiScaleKernel<double>();

  testNearestNeighbourDistance<float> ();
  testNearestNeighbourDistance<double>();

  testWassersteinDistance<float> ();
  testWassersteinDistance<double>();
}
