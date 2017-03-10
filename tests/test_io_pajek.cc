#include "config/Base.hh"

#include "tests/Base.hh"

#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include "topology/io/Pajek.hh"

#include <set>

template <class D, class V> void test( const std::string& filename )
{
  ALEPH_TEST_BEGIN( "Pajek file parsing" );

  using Simplex           = aleph::topology::Simplex<D, V>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K;

  aleph::topology::io::PajekReader reader;
  reader( filename, K );

  ALEPH_ASSERT_EQUAL( K.size(), 22 );
  ALEPH_ASSERT_EQUAL( std::count_if( K.begin(), K.end(), [] (const Simplex& s) { return s.dimension() == 0; } ), 10 );
  ALEPH_ASSERT_EQUAL( std::count_if( K.begin(), K.end(), [] (const Simplex& s) { return s.dimension() == 1; } ), 12 );

  std::set<V> vertices;
  K.vertices( std::inserter( vertices, vertices.begin() ) );

  ALEPH_ASSERT_EQUAL( *vertices.begin(),  1  );
  ALEPH_ASSERT_EQUAL( *vertices.rbegin(), 10 );

  ALEPH_TEST_END();
}

int main()
{
  std::vector<std::string> inputs = {
    CMAKE_SOURCE_DIR + std::string( "/tests/input/Simple.net" ),
    CMAKE_SOURCE_DIR + std::string( "/tests/input/Simple_with_labels.net" )
  };

  for( auto&& input : inputs )
  {
    test<double,unsigned>      ( input );
    test<double,unsigned short>( input );
    test<float, unsigned>      ( input );
    test<float, unsigned short>( input );
  }
}
