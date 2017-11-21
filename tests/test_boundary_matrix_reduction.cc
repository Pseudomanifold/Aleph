#include <tests/Base.hh>

#include <aleph/persistentHomology/Calculation.hh>
#include <aleph/persistentHomology/algorithms/Standard.hh>
#include <aleph/persistentHomology/algorithms/Twist.hh>

#include <aleph/topology/BoundaryMatrix.hh>

#include <aleph/topology/representations/Set.hh>
#include <aleph/topology/representations/Vector.hh>

#include <vector>

template <class T> void testNonSquare()
{
  ALEPH_TEST_BEGIN( "Boundary matrix reduction for non-square matrices" );

  using namespace aleph;
  using namespace topology;
  using namespace representations;

  using Representation = Vector<T>;
  using Matrix         = BoundaryMatrix<Representation>;
  using Index          = typename Matrix::Index;

  // Rectangular matrix reduction --------------------------------------

  Matrix M;
  M.setNumColumns( Index(4) );

  std::vector<T> columnA = { 1+4, 3+4, 4+4      };
  std::vector<T> columnB = { 2+4, 7+4, 8+4      };
  std::vector<T> columnC = { 5+4, 4+4, 7+4      };
  std::vector<T> columnE = { 1+4, 2+4, 5+4, 6+4 };

  M.setColumn( 0, columnA.begin(), columnA.end() );
  M.setColumn( 1, columnB.begin(), columnB.end() );
  M.setColumn( 2, columnC.begin(), columnC.end() );
  M.setColumn( 3, columnE.begin(), columnE.end() );

  using StandardAlgorithm = aleph::persistentHomology::algorithms::StandardRectangular;
  StandardAlgorithm algorithm;

  algorithm( M );

  Index i;
  bool valid;

  std::tie( i, valid ) = M.getMaximumIndex( 3 );

  ALEPH_ASSERT_EQUAL( valid, true  );
  ALEPH_ASSERT_EQUAL( i    ,   6+4 );

  // Quadratic matrix reduction w/ constraints -------------------------

  Matrix N;
  N.setNumColumns( Index(12) );

  columnA = { 1+3, 3+3, 4+3      };
  columnB = { 2+3, 7+3, 8+3      };
  columnC = { 5+3, 4+3, 7+3      };
  columnE = { 1+3, 2+3, 5+3, 6+3 };

  N.setColumn( 0, columnA.begin(), columnA.end() );
  N.setColumn( 1, columnB.begin(), columnB.end() );
  N.setColumn( 2, columnC.begin(), columnC.end() );
  N.setColumn( 3, columnE.begin(), columnE.end() );

  auto N1 = N;

  algorithm( N1 );

  std::tie( i, valid ) = N1.getMaximumIndex( 3 );

  ALEPH_ASSERT_EQUAL( valid, true  );
  ALEPH_ASSERT_EQUAL( i    ,   6+3 );

  {
    auto pairing = calculatePersistencePairing( N, false );
    std::cerr << "PAIRING: " << pairing << "\n";
  }

  std::cerr << std::string( 72, '-' ) << "\n"
            << "Original space:\n"
            << std::string( 72, '-' ) << "\n\n";

  for( Index j = 0; j < N1.getNumColumns(); j++ )
  {
    std::tie( i, valid ) = N1.getMaximumIndex( j );

    if( valid )
      std::cerr << j << ": " << i << "\n";
    else
      std::cerr << j << ": -" << "\n";
  }


  auto N2 = N.dualize();

  algorithm( N2 );

  std::cerr << N2 << "\n";

  std::cerr << std::string( 72, '-' ) << "\n"
            << "Dual space:\n"
            << std::string( 72, '-' ) << "\n\n";

  for( Index j = 0; j < N2.getNumColumns(); j++ )
  {
    std::tie( i, valid ) = N2.getMaximumIndex( j );

    if( valid )
      std::cerr << 12 - 1 - i << ": " << 12 - 1 - j << "\n";
    else
      std::cerr << 12 - 1 - j << ": -" << "\n";
  }

  auto N3      = N.dualize();
  auto pairing = calculatePersistencePairing( N3, false );

  // Quadratic matrix with re-ordering ---------------------------------

  Matrix O;
  O.setNumColumns( Index(12) );

  columnA = { 1-1, 5-1, 6-1 };
  columnB = { 2-1, 7-1, 8-1 };
  columnC = { 3-1, 6-1, 7-1 };
  columnE = { 4-1, 5-1, 8-1 };

  O.setColumn( 8, columnA.begin(), columnA.end() );
  O.setColumn( 9, columnB.begin(), columnB.end() );
  O.setColumn(10, columnC.begin(), columnC.end() );
  O.setColumn(11, columnE.begin(), columnE.end() );

  pairing = calculatePersistencePairing( O, false );

  std::cerr << "Pairing (quadratic, with re-ordering):\n";

  for( auto&& pair : pairing )
    if( pair.first <= 4 )
      std::cerr << pair.first << ": " << pair.second << "\n";

  ALEPH_TEST_END();
}

template <class M> void reduceBoundaryMatrix( const M& m )
{
  ALEPH_TEST_BEGIN( "Boundary matrix reduction" );

  ALEPH_ASSERT_THROW( m.getNumColumns() > 0 );

  using StandardAlgorithm = aleph::persistentHomology::algorithms::Standard;
  using TwistAlgorithm    = aleph::persistentHomology::algorithms::Twist;

  using Index   = typename M::Index;
  using Pairing = aleph::PersistencePairing<Index>;

  std::vector<Pairing> pairings;
  pairings.reserve( 4 );

  pairings.push_back( aleph::calculatePersistencePairing<StandardAlgorithm>( m ) );
  pairings.push_back( aleph::calculatePersistencePairing<StandardAlgorithm>( m.dualize() ) );

  pairings.push_back( aleph::calculatePersistencePairing<TwistAlgorithm>( m ) );
  pairings.push_back( aleph::calculatePersistencePairing<TwistAlgorithm>( m.dualize() ) );

  ALEPH_ASSERT_THROW( m != m.dualize() );
  ALEPH_ASSERT_THROW( m == m.dualize().dualize() );

  for( auto&& pairing : pairings )
  {
    ALEPH_ASSERT_THROW( pairing.empty() == false );
    ALEPH_ASSERT_THROW( pairing.size()  == 4 );
  }

  for( auto&& pairing1 : pairings )
    for( auto&& pairing2 : pairings )
      ALEPH_ASSERT_THROW( pairing1 == pairing2 );

  ALEPH_ASSERT_THROW( pairings.front().contains( Index(0) ) );
  ALEPH_ASSERT_THROW( pairings.front().contains( Index(1), Index(3) ) );
  ALEPH_ASSERT_THROW( pairings.front().contains( Index(2), Index(4) ) );
  ALEPH_ASSERT_THROW( pairings.front().contains( Index(5), Index(6) ) );

  ALEPH_TEST_END();
}

template <class T> void setupBoundaryMatrix()
{
  using namespace aleph;
  using namespace topology;
  using namespace representations;

  ALEPH_TEST_BEGIN( "Boundary matrix setup & loading" );

  using Set    = Set<T>;
  using Vector = Vector<T>;

  auto m1 = BoundaryMatrix<Set>::load( CMAKE_SOURCE_DIR + std::string( "/tests/input/Triangle.txt" ) );
  auto m2 = BoundaryMatrix<Vector>::load( CMAKE_SOURCE_DIR + std::string( "/tests/input/Triangle.txt" ) );

  reduceBoundaryMatrix( m1 );
  reduceBoundaryMatrix( m2 );

  ALEPH_TEST_END();
}

int main()
{
  setupBoundaryMatrix<unsigned int> ();
  setupBoundaryMatrix<unsigned long>();
  setupBoundaryMatrix<int>();
  setupBoundaryMatrix<long>();

  testNonSquare<int>          ();
  testNonSquare<long>         ();
  testNonSquare<unsigned int> ();
  testNonSquare<unsigned long>();
}
