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

  using StandardAlgorithm = aleph::persistentHomology::algorithms::Standard;
  StandardAlgorithm algorithm;
  algorithm( M );

  std::cout << M << "\n";

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
