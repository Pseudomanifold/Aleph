#include <tests/Base.hh>

#include <aleph/persistenceDiagrams/Mean.hh>
#include <aleph/persistenceDiagrams/Norms.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/PersistenceIndicatorFunction.hh>

#include <aleph/persistenceDiagrams/distances/Bottleneck.hh>
#include <aleph/persistenceDiagrams/distances/Hausdorff.hh>
#include <aleph/persistenceDiagrams/distances/NearestNeighbour.hh>
#include <aleph/persistenceDiagrams/distances/Wasserstein.hh>

#include <aleph/persistenceDiagrams/kernels/KernelEmbedding.hh>
#include <aleph/persistenceDiagrams/kernels/MultiScaleKernel.hh>

#include <algorithm>
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

template <class T> void testBottleneckDistance()
{
  ALEPH_TEST_BEGIN( "Bottleneck distance" );

  using Diagram = aleph::PersistenceDiagram<T>;
  using namespace aleph::distances;

  Diagram D1;
  D1.add( T(0.9), T(1.0) );
  D1.add( T(1.9), T(2.0) );
  D1.add( T(2.9), T(3.0) );
  D1.add( T(3.9), T(4.0) );

  {
    auto d11 = bottleneckDistance( D1, D1 );

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
    auto d12 = bottleneckDistance( D1, D2 );
    auto d21 = bottleneckDistance( D2, D1 );

    ALEPH_ASSERT_THROW( d12 > T() );
    ALEPH_ASSERT_THROW( d21 > T() );

    ALEPH_ASSERT_EQUAL( d12, d21 );
    ALEPH_ASSERT_EQUAL( d21, T(9.9)-T(4.0) );
  }

  ALEPH_TEST_END();
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

template <class T> void testHausdorffDistance()
{
  using PersistenceDiagram = aleph::PersistenceDiagram<T>;

  ALEPH_TEST_BEGIN( "Hausdorff distance" );

  auto pd = createRandomPersistenceDiagram<T>( 25 );
  auto d0 = aleph::distances::hausdorffDistance( pd, pd );
  auto d1 = aleph::distances::hausdorffDistance( PersistenceDiagram(), PersistenceDiagram() );
  auto d2 = aleph::distances::hausdorffDistance( pd                  , PersistenceDiagram() );
  auto d3 = aleph::distances::hausdorffDistance( PersistenceDiagram(), pd );

  ALEPH_ASSERT_EQUAL( d0, T() );
  ALEPH_ASSERT_EQUAL( d1, T() );
  ALEPH_ASSERT_EQUAL( d2, d3 );
  ALEPH_ASSERT_EQUAL( d2, std::numeric_limits<T>::infinity() );

  ALEPH_TEST_END();
}

template <class T> void testPersistenceIndicatorFunction()
{
  ALEPH_TEST_BEGIN( "Persistence indicator function" );

  using PersistenceDiagram = aleph::PersistenceDiagram<T>;
  using StepFunction       = aleph::math::StepFunction<T>;

  unsigned numSamples = 20;
  unsigned sampleSize = 50;

  std::vector<PersistenceDiagram> diagrams;
  diagrams.reserve( numSamples );

  for( unsigned i = 0; i < numSamples; i++ )
    diagrams.emplace_back( createRandomPersistenceDiagram<T>( sampleSize ) );

  std::vector<StepFunction> indicatorFunctions;
  indicatorFunctions.reserve( numSamples );

  for( unsigned i = 0; i < numSamples; i++ )
    indicatorFunctions.emplace_back( aleph::persistenceIndicatorFunction( diagrams.at(i) ) );

  // FIXME: not sure what to do with these values...
  for( unsigned i = 0; i < numSamples; i++ )
  {
    auto&& D1 = diagrams.at(i);
    auto&& f1 = indicatorFunctions.at(i);

    for( unsigned j = i+1; j < numSamples; j++ )
    {
      auto&& D2 = diagrams.at(j);
      auto&& f2 = indicatorFunctions.at(j);

      auto h  = aleph::distances::hausdorffDistance( D1, D2 );
      auto w1 = aleph::distances::wassersteinDistance( D1, D2, T(1) );
      auto d1 = (f1-f2).abs().integral();
      auto p1 = aleph::totalPersistence( D1, 1.0 );
      auto p2 = aleph::totalPersistence( D2, 1.0 );

      std::cout << h << "," << w1 << "," << d1 << "," << p1 << "," << p2 << "\n";
    }
  }

  // Transform persistence diagrams to prove that the bound given by
  // the maximum of the two persistence values is strict.
  for( unsigned i = 0; i < numSamples; i++ )
  {
    auto&& D1 = diagrams.at(i);
    auto   D2 = D1;

    T offset = T(10);
    std::transform( D2.begin(), D2.end(), D2.begin(),
                    [&offset] ( const typename PersistenceDiagram::Point& p )
                    {
                      return typename PersistenceDiagram::Point( p.x() + offset, p.y() + offset );
                    }
    );

    auto p1 = aleph::totalPersistence( D1, 1 );
    auto p2 = aleph::totalPersistence( D2, 1 );

    ALEPH_ASSERT_THROW( std::abs( p1 - p2 ) < 1e-4 );

    auto f = aleph::persistenceIndicatorFunction( D1 );
    auto g = aleph::persistenceIndicatorFunction( D2 );

    auto d1 = (f-g).abs().integral();
    auto d2 = (g-f).abs().integral();

    ALEPH_ASSERT_THROW( d1 >= p1+p2 );
    ALEPH_ASSERT_THROW( d2 >= p1+p2 );
  }

  ALEPH_TEST_END();
}

template <class T> void testKernelEmbedding()
{
  ALEPH_TEST_BEGIN( "Kernel embedding" );

  auto D1 = createRandomPersistenceDiagram<T>( 50 );
  auto D2 = createRandomPersistenceDiagram<T>( 50 );

  auto d1 = aleph::gaussianKernel( D1, D2, 1.0, 1.0, 1.0, 1.0 );
  auto d2 = aleph::gaussianKernel( D1, D2, 1.0, 1.0, 2.0, 2.0 );

  // Non-negativity
  ALEPH_ASSERT_THROW( d1 > 0.0 );
  ALEPH_ASSERT_THROW( d2 > 0.0 );

  // Symmetry
  ALEPH_ASSERT_EQUAL( aleph::gaussianKernel(D1, D1, 1.0, 1.0, 1.0, 1.0), aleph::gaussianKernel(D2, D2, 1.0, 1.0, 1.0, 1.0) );
  ALEPH_ASSERT_EQUAL( aleph::gaussianKernel(D1, D1, 1.0, 1.0, 1.0, 1.0), 1.0 );

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
  ALEPH_TEST_BEGIN( "Wasserstein distance" );

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

    ALEPH_ASSERT_EQUAL( d12, d21 );
    ALEPH_ASSERT_THROW( std::abs( d12 -  T( 3.05 ) ) < 1e-8 );
  }

  ALEPH_TEST_END();
}

int main(int, char**)
{
  testBottleneckDistance<float> ();
  testBottleneckDistance<double>();

  testFrechetMean<float> ();
  testFrechetMean<double>();

  testHausdorffDistance<float> ();
  testHausdorffDistance<double>();

  testKernelEmbedding<float> ();
  testKernelEmbedding<double>();

  testMultiScaleKernel<float> ();
  testMultiScaleKernel<double>();

  testNearestNeighbourDistance<float> ();
  testNearestNeighbourDistance<double>();

  testPersistenceIndicatorFunction<float> ();
  testPersistenceIndicatorFunction<double>();

  testWassersteinDistance<float> ();
  testWassersteinDistance<double>();
}
