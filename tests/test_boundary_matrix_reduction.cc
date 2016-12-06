#include "algorithms/Standard.hh"
#include "algorithms/Twist.hh"

#include "config/Base.hh"

#include "tests/Base.hh"

#include "persistentHomology/Calculation.hh"

#include "topology/BoundaryMatrix.hh"

#include "representations/Set.hh"
#include "representations/Vector.hh"

#include <vector>

template <class M> void reduceBoundaryMatrix( const M& m )
{
  ALEPH_TEST_BEGIN( "Boundary matrix reduction" );

  ALEPH_ASSERT_THROW( m.getNumColumns() > 0 );

  using StandardAlgorithm = aleph::algorithms::Standard;
  using TwistAlgorithm    = aleph::algorithms::Twist;

  using Index   = typename M::Index;
  using Pairing = aleph::PersistencePairing<Index>;

  std::vector<Pairing> pairings;
  pairings.reserve( 4 );

  pairings.push_back( aleph::calculatePersistencePairing<StandardAlgorithm>( m ) );
  pairings.push_back( aleph::calculatePersistencePairing<StandardAlgorithm>( m.dualize() ) );

  pairings.push_back( aleph::calculatePersistencePairing<TwistAlgorithm>( m ) );
  pairings.push_back( aleph::calculatePersistencePairing<TwistAlgorithm>( m.dualize() ) );

  for( auto&& pairing : pairings )
  {
    ALEPH_ASSERT_THROW( pairing.empty() == false );
    ALEPH_ASSERT_THROW( pairing.size()  == 3 );
  }

  for( auto&& pairing1 : pairings )
    for( auto&& pairing2 : pairings )
      ALEPH_ASSERT_THROW( pairing1 == pairing2 );

  ALEPH_ASSERT_THROW( pairings.front().contains( Index(1), Index(3) ) );
  ALEPH_ASSERT_THROW( pairings.front().contains( Index(2), Index(4) ) );
  ALEPH_ASSERT_THROW( pairings.front().contains( Index(5), Index(6) ) );

  ALEPH_TEST_END();
}

template <class T> void setupBoundaryMatrix()
{
  ALEPH_TEST_BEGIN( "Boundary matrix setup & loading" );

  using Set    = aleph::representations::Set<T>;
  using Vector = aleph::representations::Vector<T>;

  auto m1 = aleph::BoundaryMatrix<Set>::load( CMAKE_SOURCE_DIR + std::string( "/tests/input/Triangle.txt" ) );
  auto m2 = aleph::BoundaryMatrix<Vector>::load( CMAKE_SOURCE_DIR + std::string( "/tests/input/Triangle.txt" ) );

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
}
