#include <tests/Base.hh>
#include <aleph/topology/RandomGraph.hh>

#include <algorithm>

using Simplex = aleph::topology::Simplex<short, unsigned>;

auto isVertex = [] ( const Simplex& s )
{
  return s.dimension() == 0;
};

void testERG()
{
  ALEPH_TEST_BEGIN( "Erdos--Renyi graph" );

  auto K0 = aleph::topology::generateErdosRenyiGraph( 10, 0.0 );
  auto K1 = aleph::topology::generateErdosRenyiGraph( 10, 0.1 );
  auto K2 = aleph::topology::generateErdosRenyiGraph( 10, 0.5 );
  auto K3 = aleph::topology::generateErdosRenyiGraph( 10, 1.0 );

  ALEPH_ASSERT_EQUAL( std::count_if( K0.begin(), K0.end(), isVertex ), 10 );
  ALEPH_ASSERT_EQUAL( std::count_if( K1.begin(), K1.end(), isVertex ), 10 );
  ALEPH_ASSERT_EQUAL( std::count_if( K2.begin(), K2.end(), isVertex ), 10 );
  ALEPH_ASSERT_EQUAL( std::count_if( K3.begin(), K3.end(), isVertex ), 10 );

  ALEPH_ASSERT_THROW( K0.size() <  K3.size() );
  ALEPH_ASSERT_THROW( K1.size() <= K3.size() );
  ALEPH_ASSERT_THROW( K2.size() <= K3.size() );

  ALEPH_ASSERT_EQUAL( K0.size(), 10 );
  ALEPH_ASSERT_EQUAL( K3.size(), 10 + 10*9 / 2 );

  ALEPH_TEST_END();
}

void testWRG()
{
  ALEPH_TEST_BEGIN( "Weighted random graph" );

  auto K0 = aleph::topology::generateWeightedRandomGraph( 10, 0.0 );
  auto K1 = aleph::topology::generateWeightedRandomGraph( 10, 0.5 );

  ALEPH_ASSERT_EQUAL( K0.size(), 10 );

  ALEPH_ASSERT_THROW( K0.size() <= K1.size() );

  if( K1.size() > 10 )
  {
    unsigned maxEdgeWeight = 0;

    for( auto&& s : K1 )
      maxEdgeWeight = std::max( maxEdgeWeight, s.data() );

    ALEPH_ASSERT_THROW( maxEdgeWeight >= 1 );
  }

  ALEPH_TEST_END();
}


int main( int, char** )
{
  testERG();
  testWRG();
}
