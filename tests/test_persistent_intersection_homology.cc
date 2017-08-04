#include <tests/Base.hh>

#include <aleph/containers/PointCloud.hh>

#include <aleph/persistentHomology/algorithms/Standard.hh>

#include <aleph/persistentHomology/Calculation.hh>
#include <aleph/persistentHomology/PhiPersistence.hh>

#include <aleph/topology/BarycentricSubdivision.hh>
#include <aleph/topology/Conversions.hh>
#include <aleph/topology/QuotientSpaces.hh>
#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <random>
#include <unordered_map>
#include <vector>

#include <cmath>

template <class T> aleph::containers::PointCloud<T> sampleFromDisk( T r, unsigned n )
{
  std::random_device rd;
  std::mt19937 rng( rd() );

  std::uniform_real_distribution<T> rDistribution  ( 0, std::next( r, std::numeric_limits<T>::max() ) );
  std::uniform_real_distribution<T> phiDistribution( 0, static_cast<T>( 2 * M_PI ) );

  aleph::containers::PointCloud<T> pc( n, 2 );

  for( unsigned i = 0; i < n; i++ )
  {
    auto phi = phiDistribution( rng );
    auto r   = rDistribution( rng );
    auto x   = r * std::cos( phi );
    auto y   = r * std::sin( phi );

    pc.set( i, {x,y} );
  }

  return pc;
}

template <class T> aleph::containers::PointCloud<T> createSpokes( T r, unsigned n, unsigned K )
{
  std::random_device rd;
  std::mt19937 rng( rd() );

  std::uniform_real_distribution<T> phiDistribution( 0, static_cast<T>( 2 * M_PI ) );
  aleph::containers::PointCloud<T> pc( n*K, 2 );

  for( unsigned i = 0; i < n; i++ )
  {
    auto phi = phiDistribution( rng );
    auto x0  = r * std::cos( phi );
    auto y0  = r * std::sin( phi );
    auto x1  = x0;
    auto y1  = y0;

    for( unsigned k = 0; k < K; k++ )
    {
      pc.set( K*i+k, {x1,y1} );

      x1 += 0.05 * x0;
      y1 += 0.05 * x0;
    }
  }

  return pc;
}

template <class T> aleph::containers::PointCloud<T> makeDiskWithFlares()
{
  auto pcDisk   = sampleFromDisk( T(1), 300      );
  auto pcFlares = createSpokes(   T(1),   3, 100 );

  ALEPH_ASSERT_EQUAL( pcDisk.dimension(), pcFlares.dimension() );

  aleph::containers::PointCloud<T> pc( pcDisk.size() + pcFlares.size(), pcDisk.dimension() );
  return pc;
}



template <class T> void test()
{
  ALEPH_TEST_BEGIN( "Persistent intersection homology: simple example" );

  using Simplex           = aleph::topology::Simplex<T>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  std::vector<Simplex> simplices =
  {
    {0},
    {1},
    {2},
    {3},
    {4},
    {0,1}, {0,3}, {0,4},
    {1,2}, {1,4},
    {2,3}, {2,4},
    {3,4},
    {0,3,4}, // A
    {1,2,4}, // B
    {2,3,4}, // C
    {0,1,4}  // E
  };

  std::map<Simplex, bool> phi;

  for( auto&& simplex : simplices )
  {
    if( simplex.contains(4) == false || simplex.dimension() == 2 )
      phi[simplex] = true;
    else
      phi[simplex] = false;
  }

  SimplicialComplex K( simplices.begin(), simplices.end() );
  SimplicialComplex L;

  std::size_t s    = 0;
  std::tie( L, s ) =
    aleph::partition( K,
                      [&phi] ( const Simplex& s )
                      {
                        return phi.at(s);
                      } );

  ALEPH_ASSERT_EQUAL( K.size(), L.size() );

  auto boundaryMatrix = aleph::topology::makeBoundaryMatrix( L, s );
  auto indexA         = L.index( {0,3,4} );

  using IndexType     = typename decltype(boundaryMatrix)::Index;
  auto columnA        = boundaryMatrix.getColumn( static_cast<IndexType>( indexA ) );

  ALEPH_ASSERT_EQUAL( columnA.size(), 3 );

  aleph::persistentHomology::algorithms::Standard algorithm;
  algorithm( boundaryMatrix );

  unsigned numAllowableChains    = 0;
  unsigned numAllowableTwoChains = 0;

  for( IndexType i = 0; i < boundaryMatrix.getNumColumns(); i++ )
  {
    auto lowestOne = boundaryMatrix.getMaximumIndex( i );
    if( lowestOne.second && lowestOne.first <= s )
    {
      ++numAllowableChains;

      if( L.at(i).dimension() == 2 )
        ++numAllowableTwoChains;
    }
  }

  ALEPH_ASSERT_THROW( numAllowableChains >= numAllowableTwoChains );
  ALEPH_ASSERT_EQUAL( numAllowableTwoChains, 1 );

  ALEPH_TEST_END();
}

template <class T> void testCircleWithWhisker()
{
  ALEPH_TEST_BEGIN( "Persistent intersection homology: circle plus whisker" );

  using Simplex           = aleph::topology::Simplex<T>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K = {
    {0},
    {1},
    {2},
    {0,1}, {0,2},
    {1,2}
  };

  SimplicialComplex L =
  {
    {0},
    {1},
    {2},
    {3},
    {0,1}, {0,2}, {0,3},
    {1,2}
  };

  K.sort();
  L.sort();

  {
    auto D1 = aleph::calculatePersistenceDiagrams( K );
    auto D2 = aleph::calculatePersistenceDiagrams( L );

    ALEPH_ASSERT_EQUAL( D1.size(),          D2.size() );
    ALEPH_ASSERT_EQUAL( D1.front().betti(), D2.front().betti() );
  }

  SimplicialComplex X0 = { {0}, {1}, {2} };
  SimplicialComplex X1 = K;

  SimplicialComplex Y0 = { {0} };
  SimplicialComplex Y1 = L;

  auto D1 = aleph::calculateIntersectionHomology( K, {X0,X1}, aleph::Perversity( {-1} ) );
  auto D2 = aleph::calculateIntersectionHomology( L, {Y0,Y1}, aleph::Perversity( {-1} ) );

  ALEPH_ASSERT_THROW( D1.empty() == false );
  ALEPH_ASSERT_THROW( D2.empty() == false );
  ALEPH_ASSERT_EQUAL( D1.front().dimension(), 0 );
  ALEPH_ASSERT_EQUAL( D2.front().dimension(), 0 );

  ALEPH_ASSERT_EQUAL( D1.front().betti(), 1 );
  ALEPH_ASSERT_EQUAL( D2.front().betti(), 2 );

  ALEPH_TEST_END();
}
template <class T> void testQuotientSpaces()
{
  ALEPH_TEST_BEGIN( "Persistent intersection homology: quotient spaces" );

  using Simplex           = aleph::topology::Simplex<T>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K = {
    {0},
    {1},
    {2},
    {3},
    {0,1}, {0,2}, {0,3},
    {1,2}, {1,3},
    {2,3},
    {0,1,2}, {0,1,3}, {0,2,3},
    {1,2,3}
  };

  auto C = aleph::topology::cone( K );
  auto S = aleph::topology::suspension( K );

  ALEPH_ASSERT_THROW( C.empty() == false );
  ALEPH_ASSERT_THROW( S.empty() == false );
  ALEPH_ASSERT_EQUAL( C.size(), 2 * K.size() + 1);
  ALEPH_ASSERT_EQUAL( S.size(), 3 * K.size() + 2);

  S.sort();

  bool dualize                    = true;
  bool includeAllUnpairedCreators = true;

  auto D1 = aleph::calculatePersistenceDiagrams(K, dualize, includeAllUnpairedCreators);
  auto D2 = aleph::calculatePersistenceDiagrams(S, dualize, includeAllUnpairedCreators);

  ALEPH_ASSERT_EQUAL( D1.size(), 3 );
  ALEPH_ASSERT_EQUAL( D2.size(), 4 );

  std::vector<std::size_t> expectedBettiNumbersK = {1,0,1};
  std::vector<std::size_t> expectedBettiNumbersS = {1,0,0,1};

  std::vector<std::size_t> bettiNumbersK;
  std::vector<std::size_t> bettiNumbersS;

  using PersistenceDiagram = decltype( D1.front() );

  std::transform( D1.begin(), D1.end(), std::back_inserter( bettiNumbersK ), [] ( const PersistenceDiagram& D ) { return D.betti(); } );
  std::transform( D2.begin(), D2.end(), std::back_inserter( bettiNumbersS ), [] ( const PersistenceDiagram& D ) { return D.betti(); } );

  ALEPH_ASSERT_THROW( bettiNumbersK == expectedBettiNumbersK );
  ALEPH_ASSERT_THROW( bettiNumbersS == expectedBettiNumbersS );

  ALEPH_TEST_END();
}

template <class T> void testSphere()
{
  ALEPH_TEST_BEGIN( "Persistent intersection homology: sphere triangulation" );

  using Simplex           = aleph::topology::Simplex<T>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K = {
    {0},
    {1},
    {2},
    {3},
    {0,1}, {0,2}, {0,3},
    {1,2}, {1,3},
    {2,3},
    {0,1,2}, {0,1,3}, {0,2,3},
    {1,2,3}
  };

  SimplicialComplex X0 = { {0}, {1}, {2}, {3} };
  SimplicialComplex X1 = K;

  auto D1 = aleph::calculateIntersectionHomology( K, {X0,X1}, aleph::Perversity( {0,0} ) );

  // This demonstrates that the triangulation does not have any
  // allowable vertices. Hence, no intersection homology exists
  // in dimension 0.
  ALEPH_ASSERT_EQUAL( D1.size(),              1 );
  ALEPH_ASSERT_EQUAL( D1.front().dimension(), 2 );

  SimplicialComplex L = K;

  {
    aleph::topology::BarycentricSubdivision subdivision;
    L = subdivision( L );
    L.sort();
  }

  auto Y1 = L;
  auto D2 = aleph::calculateIntersectionHomology( L, {X0,K}, aleph::Perversity( {0,0} ) );

  // This demonstrates that the barycentric subdivision of the space,
  // i.e. another triangulation, may influence the results.
  ALEPH_ASSERT_EQUAL( D2.size(),              3 );
  ALEPH_ASSERT_EQUAL( D2.front().dimension(), 0 );
  ALEPH_ASSERT_EQUAL( D2.front().betti(),     1 );

  ALEPH_TEST_END();
}

template <class T> void testWedgeOfTwoCircles()
{
  ALEPH_TEST_BEGIN( "Persistent intersection homology: wedge of two circles" );

  using Simplex           = aleph::topology::Simplex<T>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K = {
    {0},
    {1},
    {2},
    {3},
    {4},
    {5},
    {6},
    {0,1}, {0,6},
    {1,2},
    {2,3}, {2,5}, {2,6},
    {3,4},
    {4,5}
  };

  SimplicialComplex X0 = { {2} };
  SimplicialComplex X1 = K;

  // This example demonstrates the dependence on the filtration or
  // rather the stratification of the complex.
  //
  // Using the equal perversity as for the previous example, a new
  // component is being created.
  SimplicialComplex Y0 = { {0}, {2} };
  SimplicialComplex Y1 = K;

  auto D1 = aleph::calculateIntersectionHomology( K, {X0,X1}, aleph::Perversity( {-1} ) );
  auto D2 = aleph::calculateIntersectionHomology( K, {X0,X1}, aleph::Perversity( { 0} ) );
  auto D3 = aleph::calculateIntersectionHomology( K, {Y0,Y1}, aleph::Perversity( {-1} ) );

  ALEPH_ASSERT_EQUAL( D1.size(), 1 );
  ALEPH_ASSERT_EQUAL( D2.size(), 2 );
  ALEPH_ASSERT_EQUAL( D3.size(), 1 );

  ALEPH_ASSERT_EQUAL( D1[0].betti(), 2 );
  ALEPH_ASSERT_EQUAL( D3[0].betti(), 3 );

  // TODO: is this correct? In his Ph.D. thesis "Analyzing Stratified
  // Spaces Using Persistent Versions of Intersection and Local
  // Homology", Bendich states that this should be 0...
  ALEPH_ASSERT_EQUAL( D2[0].betti(), 1 );
  ALEPH_ASSERT_EQUAL( D2[1].betti(), 2 );

  ALEPH_TEST_END();
}

int main(int, char**)
{
  test<float> ();
  test<double>();

  testCircleWithWhisker<float> ();
  testCircleWithWhisker<double>();

  testQuotientSpaces<float> ();
  testQuotientSpaces<double>();

  testSphere<float> ();
  testSphere<double>();

  testWedgeOfTwoCircles<float> ();
  testWedgeOfTwoCircles<double>();
}
