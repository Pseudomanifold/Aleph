#include <tests/Base.hh>

#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/CechComplex.hh>
#include <aleph/geometry/VietorisRipsComplex.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <aleph/persistentHomology/Calculation.hh>
#include <aleph/persistentHomology/PhiPersistence.hh>

#include <aleph/topology/BarycentricSubdivision.hh>
#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>
#include <aleph/topology/Skeleton.hh>
#include <aleph/topology/Spine.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/topology/io/LinesAndPoints.hh>

#include <random>
#include <vector>

#include <cmath>

template <class T> void testDisk()
{
  ALEPH_TEST_BEGIN( "Spine: disk" );

  using DataType   = bool;
  using VertexType = T;

  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  std::vector<Simplex> simplices;

  unsigned n = 7;
  for( unsigned i = 0; i < n; i++ )
  {
    if( i+1 < n )
      simplices.push_back( Simplex( {T(0),T(i+1),T(i+2)} ) );
    else
      simplices.push_back( Simplex( {T(0),T(i+1),T(  1)} ) );
  }

  SimplicialComplex K( simplices.begin(), simplices.end() );

  K.createMissingFaces();
  K.sort();

  auto L = aleph::topology::spine( K );
  auto M = aleph::topology::dumb::spine( K );

  ALEPH_ASSERT_THROW( L.size() < K.size() );
  ALEPH_ASSERT_EQUAL( L.size(), 1 );
  ALEPH_ASSERT_EQUAL( M.size(), 1 );

  // Note that it does not make sense to check whether both spines
  // resulted in the *same* vertex. Every vertex is equally likely
  // to be left over; and every result is equally valid.

  ALEPH_TEST_END();
}

template <class T> void testPinchedTorus()
{
  using DataType   = T;
  using PointCloud = aleph::containers::PointCloud<DataType>;

  ALEPH_TEST_BEGIN( "Spine: pinched torus" );

  unsigned n =  40;
  unsigned m =  20;
  PointCloud pc( n*m, 3 );

  auto g = [] ( T x, T y )
  {
    return T(2) + std::sin(x/2) * std::cos(y);
  };

  std::random_device rd;
  std::mt19937 rng( rd() );

  std::normal_distribution<DataType> noise( T(0), T(0.05) );

  unsigned k = 0;
  for( unsigned i = 0; i < n; i++ )
  {
    auto x = T( 2*M_PI / n * i );
    for( unsigned j = 0; j < m; j++ )
    {
      auto y = T( 2*M_PI / m * j );

      auto x0 = g(x,y) * std::cos(x)        + noise( rng );
      auto x1 = g(x,y) * std::sin(x)        + noise( rng );
      auto x2 = std::sin(x/2) * std::sin(y) + noise( rng );

      pc.set(k++, {x0,x1,x2} );
    }
  }

  using Distance          = aleph::geometry::distances::Euclidean<DataType>;
  using NearestNeighbours = aleph::geometry::BruteForce<PointCloud, Distance>;

  auto K
    = aleph::geometry::buildVietorisRipsComplex(
      NearestNeighbours( pc ),
      DataType( 0.700 ),
      2
  );

  std::ofstream out( "/tmp/Pinched_torus.txt" );
  aleph::topology::io::LinesAndPoints lap;
  lap( out, K, pc );

  auto D1 = aleph::calculatePersistenceDiagrams( K );

  ALEPH_ASSERT_EQUAL( D1.size(), 2 );
  ALEPH_ASSERT_EQUAL( D1[0].dimension(), 0 );
  ALEPH_ASSERT_EQUAL( D1[1].dimension(), 1 );
  ALEPH_ASSERT_EQUAL( D1[1].betti(),     1 );

#if 0
  // FIXME: this is still too large to be easily processed by the
  // algorithm...

  auto L  = aleph::topology::spine( K );

  ALEPH_ASSERT_THROW( L.size() < K.size() );

  auto K0 = aleph::topology::Skeleton()(0, K);
  auto K1 = K0;
  auto K2 = K;
  auto D2 = aleph::calculateIntersectionHomology( L, {K0,K1,K2}, aleph::PerversityGM( {0} ) );

  ALEPH_ASSERT_EQUAL( D2.size(), 3 );

#endif

  ALEPH_TEST_END();
}

template <class T> void testS1vS1()
{
  using DataType   = T;
  using PointCloud = aleph::containers::PointCloud<DataType>;

  ALEPH_TEST_BEGIN( "Spine: S^1 v S^1" );

  unsigned n = 10;

  PointCloud pc( 2*n - 1, 2 );

  unsigned k = 0;
  for( unsigned i = 0; i < n; i++ )
  {
    auto x0 = DataType( std::cos( 2*M_PI / n * i ) );
    auto y0 = DataType( std::sin( 2*M_PI / n * i ) );

    if( x0 > -1 )
    {
      auto x1 = x0 + 2;
      auto y1 = y0;

      pc.set(k++, {x0, y0});
      pc.set(k++, {x1, y1});
    }

    // prevent duplication of singular point
    else
      pc.set(k++, {x0, y0} );
  }

  auto K
    = aleph::geometry::buildCechComplex(
      pc,
      DataType( 0.75 )
  );

  {
    std::ofstream out( "/tmp/SimplicialComplex.txt" );
    out << K << "\n";
  }

  {
    std::ofstream out( "/tmp/K.txt" );
    aleph::topology::io::LinesAndPoints lap;
    lap( out, K, pc );
  }

  auto D1 = aleph::calculatePersistenceDiagrams( K );

  // Persistent homology -----------------------------------------------
  //
  // This should not be surprising: it is possible to extract the two
  // circles from the data set. They form one connected component.

  ALEPH_ASSERT_THROW( D1.size() >=   2 );
  ALEPH_ASSERT_EQUAL( D1[0].betti(), 1 );
  ALEPH_ASSERT_EQUAL( D1[1].betti(), 2 );

  // Persistent intersection homology ----------------------------------
  //
  // Regardless of the stratification, it is impossible to detect the
  // singularity in dimension 0.

  auto L  = aleph::topology::BarycentricSubdivision()( aleph::topology::Skeleton()( 2, K), [] ( std::size_t dimension ) { return dimension == 0 ? 0 : 0.5; } );
  auto K0 = aleph::topology::Skeleton()( 0, K );
  auto D2 = aleph::calculateIntersectionHomology( L, {K0, aleph::topology::Skeleton()(2, K) }, aleph::Perversity( {-1} ) );

  ALEPH_ASSERT_THROW( D2.size() >=       1 );
  ALEPH_ASSERT_EQUAL( D2[0].dimension(), 0 );
  ALEPH_ASSERT_EQUAL( D2[0].betti(),     1 );

  // Spine calculation -------------------------------------------------

  auto M = aleph::topology::spine( K );
  K.sort( aleph::topology::filtrations::Data<typename decltype(K)::ValueType>() );

  {
    std::ofstream out( "/tmp/M.txt" );
    aleph::topology::io::LinesAndPoints lap;
    lap.addVertexLabels( true );
    lap( out, M, pc );
  }

  {
    std::ofstream out( "/tmp/Spine_complex.txt" );
    out << M << "\n";
  }

  {
    bool dualize                    = true;
    bool includeAllUnpairedCreators = true;

    auto D
      = aleph::calculatePersistenceDiagrams(
          M,
          dualize,
          includeAllUnpairedCreators
    );

    ALEPH_ASSERT_THROW( D.size() >=       2 );
    ALEPH_ASSERT_EQUAL( D[0].dimension(), 0 );
    ALEPH_ASSERT_EQUAL( D[1].dimension(), 1 );
    ALEPH_ASSERT_EQUAL( D[0].betti(),     1 );
    ALEPH_ASSERT_EQUAL( D[1].betti(),     2 );
  }

  ALEPH_ASSERT_THROW( M.size() < K .size() );
  ALEPH_TEST_END();

  L = aleph::topology::BarycentricSubdivision()( M, [] ( std::size_t dimension ) { return dimension == 0 ? 0 : 0.5; } );
  L.sort( aleph::topology::filtrations::Data<typename decltype(L)::ValueType>() );

  K0      = { {9}, {17} };
  auto D3 = aleph::calculateIntersectionHomology( L, {K0,M}, aleph::Perversity( {-1,0} ) );

  ALEPH_ASSERT_THROW( D3.empty() == false );
  ALEPH_ASSERT_EQUAL( D3[0].dimension(), 0 );
  ALEPH_ASSERT_EQUAL( D3[0].betti(),     3 );
}

template <class T> void testTriangle()
{
  ALEPH_TEST_BEGIN( "Triangle" );

  using DataType   = bool;
  using VertexType = T;

  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K = {
    {0,1,2},
    {0,1}, {0,2}, {1,2},
    {0}, {1}, {2}
  };

  auto L = aleph::topology::spine( K );
  auto M = aleph::topology::dumb::spine( K );

  ALEPH_ASSERT_THROW( L.size() < K.size() );
  ALEPH_ASSERT_EQUAL( L.size(), 1 );
  ALEPH_ASSERT_EQUAL( M.size(), 1 );

  ALEPH_TEST_END();
}

int main( int, char** )
{
  testDisk<short>   ();
  testDisk<unsigned>();

  testPinchedTorus<float> ();
  testPinchedTorus<double>();

  testS1vS1<float> ();
  testS1vS1<double>();

  testTriangle<short>   ();
  testTriangle<unsigned>();
}
