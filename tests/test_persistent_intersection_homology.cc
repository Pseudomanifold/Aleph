#include <tests/Base.hh>

#include <aleph/containers/DimensionalityEstimators.hh>
#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/VietorisRipsComplex.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <aleph/persistentHomology/algorithms/Standard.hh>

#include <aleph/persistentHomology/Calculation.hh>
#include <aleph/persistentHomology/PhiPersistence.hh>

#include <aleph/topology/BarycentricSubdivision.hh>
#include <aleph/topology/Conversions.hh>
#include <aleph/topology/QuotientSpaces.hh>
#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>
#include <aleph/topology/Skeleton.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <algorithm>
#include <random>
#include <unordered_map>
#include <vector>

#include <cmath>

template <class T> aleph::containers::PointCloud<T> sampleFromDisk( T r, unsigned n )
{
  std::random_device rd;
  std::mt19937 rng( rd() );

  std::uniform_real_distribution<T> rDistribution  ( 0, std::nextafter( r, std::numeric_limits<T>::max() ) );
  std::uniform_real_distribution<T> phiDistribution( 0, static_cast<T>( 2 * M_PI ) );

  aleph::containers::PointCloud<T> pc( n, 2 );

  for( unsigned i = 0; i < n; i++ )
  {
    auto phi = phiDistribution( rng );
    auto r   = std::sqrt( rDistribution( rng ) );
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

      x1 += T(0.05) * x0;
      y1 += T(0.05) * y1;
    }
  }

  return pc;
}

template <class T, class OutputIterator> aleph::containers::PointCloud<T> makeDiskWithFlares( OutputIterator result )
{
  auto pcDisk   = sampleFromDisk( T(1), 300      );
  auto pcFlares = createSpokes(   T(1),   3, 10 );

  *result++ = 300;
  *result++ = 301;
  *result++ = 310;
  *result++ = 311;
  *result++ = 320;
  *result++ = 321;

  ALEPH_ASSERT_EQUAL( pcDisk.dimension(), pcFlares.dimension() );

  return pcDisk + pcFlares;
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

  // The simplest way to model a circle using a simplicial complex, i.e.
  // the edges and vertices of a triangle.
  SimplicialComplex K = {
    {0},
    {1},
    {2},
    {0,1}, {0,2},
    {1,2}
  };

  // An additional vertex with a small 'whisker' has been added here in
  // order to show the difference between ordinary homology and
  // intersection homology.
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

  SimplicialComplex X0 = { {0} };
  SimplicialComplex X1 = K;

  SimplicialComplex Y0 = { {0} };
  SimplicialComplex Y1 = L;

  auto D1 = aleph::calculateIntersectionHomology( K, {X0,X1}, aleph::Perversity( {-1} ) );
  auto D2 = aleph::calculateIntersectionHomology( K, {X0,X1}, aleph::Perversity( { 0} ) );
  auto D3 = aleph::calculateIntersectionHomology( L, {Y0,Y1}, aleph::Perversity( {-1} ) );
  auto D4 = aleph::calculateIntersectionHomology( L, {Y0,Y1}, aleph::Perversity( { 0} ) );

  ALEPH_ASSERT_THROW( D1.empty() == false );
  ALEPH_ASSERT_THROW( D2.empty() == false );
  ALEPH_ASSERT_THROW( D3.empty() == false );
  ALEPH_ASSERT_THROW( D4.empty() == false );

  ALEPH_ASSERT_EQUAL( D1.front().dimension(), 0 );
  ALEPH_ASSERT_EQUAL( D2.front().dimension(), 0 );
  ALEPH_ASSERT_EQUAL( D3.front().dimension(), 0 );
  ALEPH_ASSERT_EQUAL( D4.front().dimension(), 0 );

  ALEPH_ASSERT_EQUAL( D1.front().betti(), 1 );
  ALEPH_ASSERT_EQUAL( D2.front().betti(), 1 );
  ALEPH_ASSERT_EQUAL( D3.front().betti(), 2 );
  ALEPH_ASSERT_EQUAL( D4.front().betti(), 1 );

  ALEPH_TEST_END();
}

template <class T> void testDiskWithFlares()
{
  ALEPH_TEST_BEGIN( "Persistent intersection homology: disk with flares" );

  using PointCloud        = aleph::containers::PointCloud<T>;
  using Distance          = aleph::geometry::distances::Euclidean<T>;
  using NearestNeighbours = aleph::geometry::BruteForce<PointCloud, Distance>;

  std::vector<std::size_t> singularIndices;

  auto pc = makeDiskWithFlares<T>( std::back_inserter( singularIndices ) );
  auto K  = aleph::geometry::buildVietorisRipsComplex(
    NearestNeighbours( pc ),
    T(0.225),
    1
  );

#if 0

  // FIXME: this is horrible...
  auto dimensionalities
    = aleph::containers::estimateLocalDimensionalityNearestNeighbours<Distance, PointCloud, NearestNeighbours>( pc,
         8,
        Distance() );

#endif

  ALEPH_ASSERT_THROW( pc.empty() == false );
  ALEPH_ASSERT_THROW( K.empty() == false );

  using SimplicialComplex = decltype(K);
  using Simplex           = typename SimplicialComplex::ValueType;
  using VertexType        = typename Simplex::VertexType;

  SimplicialComplex X0;
  for( auto&& singularIndex : singularIndices )
    X0.push_back( Simplex( static_cast<VertexType>( singularIndex ) ) );

  SimplicialComplex X1 = K;

  {
    aleph::topology::BarycentricSubdivision subdivision;
    K = subdivision( K );
    K.sort( aleph::topology::filtrations::Data<Simplex>() );
  }

  auto diagramsPH   = aleph::calculatePersistenceDiagrams( K );
  auto diagramsIH_1 = aleph::calculateIntersectionHomology( K, {X0,X1}, aleph::Perversity( {-1} ) );
  auto diagramsIH_2 = aleph::calculateIntersectionHomology( K, {X0,X1}, aleph::Perversity( { 0} ) );

  ALEPH_ASSERT_THROW( diagramsPH.empty() == false );
  ALEPH_ASSERT_EQUAL( diagramsPH.front().dimension(), 0 );

  if( diagramsPH.front().betti() != 1 )
  {
    ALEPH_TEST_END();
    return;
  }

  ALEPH_ASSERT_EQUAL( diagramsPH.front().betti()    , 1 );

  ALEPH_ASSERT_THROW( diagramsIH_1.empty() == false );
  ALEPH_ASSERT_THROW( diagramsIH_2.empty() == false );
  ALEPH_ASSERT_EQUAL( diagramsIH_1.front().dimension(), 0 );
  ALEPH_ASSERT_EQUAL( diagramsIH_2.front().dimension(), 0 );

  std::cerr << "Betti numbers:\n"
            << "  - PH  : " << diagramsPH.front().betti() << "\n"
            << "  - IH_1: " << diagramsIH_1.front().betti() << "\n"
            << "  - IH_2: " << diagramsIH_2.front().betti() << "\n";

  ALEPH_TEST_END();
}

template <class T> void testPinchedTorus()
{
  ALEPH_TEST_BEGIN( "Persistent intersection homology: pinched torus" );

  using Simplex           = aleph::topology::Simplex<T>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K = {
    {1,4,5},
    {1,2,5},
    {0,3,5},
    {0,2,5},
    {0,1,4},
    {0,3,4}
  };

  K.createMissingFaces();
  K.sort();

  using VertexType = typename Simplex::VertexType;
  K.push_back( {6} );

  {
    std::vector<Simplex> simplices;
    for( auto&& simplex : K )
    {
      if( simplex.dimension() == 1 )
      {
        auto u = simplex[0];
        auto v = simplex[1];

        // check that the simplex is only spanning a single triangle at
        // either the end or the beginning of the space
        if( std::max(u,v) <= 2 || std::min(u,v) >= 3 )
          simplices.push_back( {u,v,VertexType(6) } );
      }
    }

    for( auto&& simplex : simplices )
      K.push_back( simplex );

    K.createMissingFaces();
    K.sort();
  }

  bool dualize                    = true;
  bool includeAllUnpairedCreators = true;

  auto D1 = aleph::calculatePersistenceDiagrams(K, dualize, includeAllUnpairedCreators);

  ALEPH_ASSERT_EQUAL( D1.size(), 3 );
  ALEPH_ASSERT_EQUAL( D1[0].betti(), 1 ); // Z
  ALEPH_ASSERT_EQUAL( D1[1].betti(), 1 ); // Z
  ALEPH_ASSERT_EQUAL( D1[2].betti(), 1 ); // Z

  SimplicialComplex L = K;

  {
    aleph::topology::BarycentricSubdivision subdivision;
    L = subdivision( L );
    L.sort();
  }

  aleph::topology::Skeleton skeleton;
  auto X0 = skeleton(0, K);
  auto X1 = X0;
  auto X2 = K;

  auto D2 = aleph::calculateIntersectionHomology( L, {X0,X1,X2}, aleph::PerversityGM( {0} ) );

  ALEPH_ASSERT_EQUAL( D2.size(), 3 );
  ALEPH_ASSERT_EQUAL( D2[0].betti(), 1 );
  ALEPH_ASSERT_EQUAL( D2[1].betti(), 0 );
  ALEPH_ASSERT_EQUAL( D2[2].betti(), 1 );

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
  // in dimension 0 (likewise for dimension 1).
  ALEPH_ASSERT_EQUAL( D1.size(), 0 );

  SimplicialComplex L = K;

  {
    aleph::topology::BarycentricSubdivision subdivision;
    L = subdivision( L );
    L.sort();
  }

  auto D2 = aleph::calculateIntersectionHomology( L, {X0,K}, aleph::Perversity( {0,0} ) );

  // This demonstrates that the barycentric subdivision of the space,
  // i.e. another triangulation, may influence the results.
  ALEPH_ASSERT_EQUAL( D2.size(),              3 );
  ALEPH_ASSERT_EQUAL( D2.front().dimension(), 0 );
  ALEPH_ASSERT_EQUAL( D2.front().betti(),     1 );

  auto D3 = aleph::calculateIntersectionHomology( L, {X0,X0,K}, aleph::PerversityGM( {0} ) );

  ALEPH_ASSERT_EQUAL( D3.size(),         3 );
  ALEPH_ASSERT_EQUAL( D3[0].dimension(), 0 );
  ALEPH_ASSERT_EQUAL( D3[0].betti(),     1 );
  ALEPH_ASSERT_EQUAL( D3[1].betti(),     0 );
  ALEPH_ASSERT_EQUAL( D3[2].betti(),     1 );

  ALEPH_TEST_END();
}

template <class T> void testTorus()
{
  ALEPH_TEST_BEGIN( "Persistent intersection homology: torus" );

  using Simplex           = aleph::topology::Simplex<T>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K = {
    {0, 1, 5}, {0, 4, 5}, { 1, 2, 6}, { 1, 5, 6}, {2,3, 7}, {2, 6, 7}, {4,5,9}, {4,8,9},
    {5, 6,10}, {5, 9,10}, { 6, 7,11}, { 6,10,11}, {8,9,13}, {8,12,13},
    {9,10,14}, {9,13,14}, {10,11,15}, {10,14,15}
  };

  ALEPH_ASSERT_EQUAL( K.size(), 18 );

  // Identify faces (top and bottom) -----------------------------------

  {
    std::vector<Simplex> simplices( K.begin(), K.end() );

    std::transform( simplices.begin(), simplices.end(), simplices.begin(),
      [] ( const Simplex& s )
      {
        using VertexType = typename Simplex::VertexType;

        std::map<VertexType, VertexType> vertexMap = {
          {12,0},
          {13,1},
          {14,2},
          {15,3},
        };

        std::vector<VertexType> vertices( s.begin(), s.end() );

        for( auto&& v : vertices )
        {
          if( vertexMap.find(v) != vertexMap.end() )
            v = vertexMap.at(v);
        }

        return Simplex( vertices.begin(), vertices.end() );
      }
    );

    K = SimplicialComplex( simplices.begin(), simplices.end() );
  }

  // Identify faces (left and right) -----------------------------------

  {
    std::vector<Simplex> simplices( K.begin(), K.end() );

    std::transform( simplices.begin(), simplices.end(), simplices.begin(),
      [] ( const Simplex& s )
      {
        using VertexType = typename Simplex::VertexType;

        std::map<VertexType, VertexType> vertexMap = {
          { 3,0},
          { 7,4},
          {11,8},
        };

        std::vector<VertexType> vertices( s.begin(), s.end() );

        for( auto&& v : vertices )
        {
          if( vertexMap.find(v) != vertexMap.end() )
            v = vertexMap.at(v);
        }

        return Simplex( vertices.begin(), vertices.end() );
      }
    );

    K = SimplicialComplex( simplices.begin(), simplices.end() );
  }

  K.createMissingFaces();
  K.sort();

  bool dualize                    = true;
  bool includeAllUnpairedCreators = true;

  // Sanity check ------------------------------------------------------
  //
  // Ensures that this is a proper triangulation of a torus, at least
  // from a homological point of view.

  {
    auto D1 = aleph::calculatePersistenceDiagrams(K, dualize, includeAllUnpairedCreators);

    ALEPH_ASSERT_EQUAL( D1.size(), 3 );
    ALEPH_ASSERT_EQUAL( D1[0].betti(), 1 ); // Z
    ALEPH_ASSERT_EQUAL( D1[1].betti(), 2 ); // Z+Z
    ALEPH_ASSERT_EQUAL( D1[2].betti(), 1 ); // Z
  }

  K = aleph::topology::suspension( K );
  K.sort();

  // Sanity check ------------------------------------------------------
  //
  // Ensures that we have the expected homology groups for the suspended
  // torus data set.

  {
    auto D1 = aleph::calculatePersistenceDiagrams(K, dualize, includeAllUnpairedCreators);

    ALEPH_ASSERT_EQUAL( D1.size(), 4 );
    ALEPH_ASSERT_EQUAL( D1[0].betti(), 1 ); // Z
    ALEPH_ASSERT_EQUAL( D1[1].betti(), 0 ); // 0
    ALEPH_ASSERT_EQUAL( D1[2].betti(), 2 ); // Z+Z
    ALEPH_ASSERT_EQUAL( D1[3].betti(), 1 ); // Z
  }

  // Intersection homology ---------------------------------------------

  SimplicialComplex L = K;

  {
    aleph::topology::BarycentricSubdivision subdivision;
    L = subdivision( L );
    L.sort();
  }

  aleph::topology::Skeleton skeleton;
  auto X0 = skeleton(0, K);
  auto X1 = skeleton(1, K);
  auto X2 = X1;
  auto X3 = K;

  auto D2 = aleph::calculateIntersectionHomology( L, {X0,X1,X2,X3}, aleph::PerversityGM( {0,0} ) );
  auto D3 = aleph::calculateIntersectionHomology( L, {X0,X1,X2,X3}, aleph::PerversityGM( {0,1} ) );

  ALEPH_ASSERT_EQUAL( D2.size(), 4 );
  ALEPH_ASSERT_EQUAL( D3.size(), 4 );

  ALEPH_ASSERT_EQUAL( D2[0].betti(), 1 ); // Z
  ALEPH_ASSERT_EQUAL( D2[1].betti(), 2 ); // Z+Z
  ALEPH_ASSERT_EQUAL( D2[2].betti(), 0 ); // 0
  ALEPH_ASSERT_EQUAL( D2[3].betti(), 1 ); // Z

  ALEPH_ASSERT_EQUAL( D3[0].betti(), 1 ); // Z
  ALEPH_ASSERT_EQUAL( D3[1].betti(), 0 ); // 0
  ALEPH_ASSERT_EQUAL( D3[2].betti(), 2 ); // Z+Z
  ALEPH_ASSERT_EQUAL( D3[3].betti(), 1 ); // Z

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

template <class T> void testWeightedTriangle()
{
  ALEPH_TEST_BEGIN( "Weighted triangle" );

  using DataType          = float;
  using VertexType        = T;
  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K = {
   {0},
   {1},
   {2},
   Simplex( {0,1},   DataType(1) ), Simplex( {0,2}, DataType(2) ),
   Simplex( {1,2},   DataType(1) ),
   Simplex( {0,1,2}, DataType(2) )
  };

  aleph::topology::BarycentricSubdivision Sd;

  auto L = Sd(K, [] ( std::size_t dimension ) { return dimension == 0 ? 0 : 0.5; } );

  {
    bool useMaximum                  = true;
    bool skipOneDimensionalSimplices = true;

    L.recalculateWeights( useMaximum, skipOneDimensionalSimplices );
    L.sort( aleph::topology::filtrations::Data<Simplex>() );
  }

  ALEPH_ASSERT_THROW( L.empty() == false );
  ALEPH_ASSERT_THROW( K.size() < L.size() );
  ALEPH_ASSERT_EQUAL( L.size(), 25 );

  auto X0 = SimplicialComplex( { {0}, {1}, {2} } );
  auto X1 = X0;
  auto X2 = K;

  auto D1 = aleph::calculateIntersectionHomology( L, {X0,X1,X2}, aleph::PerversityGM( {0} ) );

  for( auto&& D : D1 )
    std::cerr << D << "\n";

  // Brute-force calculation of admissibility condition ----------------
  //
  // We have p(k) = 0 since we are using a Goresky--MacPherson
  // perversity. Moreover, only k=2 is relevant because K only
  // has 2-dimensional top simplices. This ensures that we use
  // X0 = K0, i.e the 0-skeleton of the complex.

  std::map<Simplex, bool> phi;

  for( auto&& s : L )
  {
    auto i            = s.dimension();
    auto intersection = aleph::topology::lastLexicographicalIntersection( X0, s );
    auto dimension    = intersection.empty() ? -1 : static_cast<long>( intersection.dimension() );
    bool admissible   = intersection.empty() ? true : static_cast<long>( dimension ) <= ( long(i) - long(2) );

    phi[s] = admissible;
  }

  unsigned numAdmissible0Simplices = 0;
  unsigned numAdmissible1Simplices = 0;
  unsigned numAdmissible2Simplices = 0;

  for( auto&& pair : phi )
  {
    if( pair.second )
    {
      if( pair.first.dimension() == 0 )
        ++numAdmissible0Simplices;
      else if( pair.first.dimension() == 1 )
        ++numAdmissible1Simplices;
      else
        ++numAdmissible2Simplices;
    }
  }

  ALEPH_ASSERT_EQUAL( numAdmissible0Simplices, 4 );
  ALEPH_ASSERT_EQUAL( numAdmissible1Simplices, 3 );
  ALEPH_ASSERT_EQUAL( numAdmissible2Simplices, 6 );

  ALEPH_TEST_END();
}

int main(int, char**)
{
  test<float> ();
  test<double>();

  testCircleWithWhisker<float> ();
  testCircleWithWhisker<double>();

  testDiskWithFlares<float> ();
  testDiskWithFlares<double>();

  testPinchedTorus<float> ();
  testPinchedTorus<double>();

  testQuotientSpaces<float> ();
  testQuotientSpaces<double>();

  testSphere<float> ();
  testSphere<double>();

  testTorus<float> ();
  testTorus<double>();

  testWedgeOfTwoCircles<float> ();
  testWedgeOfTwoCircles<double>();

  testWeightedTriangle<float> ();
  testWeightedTriangle<double>();
}
