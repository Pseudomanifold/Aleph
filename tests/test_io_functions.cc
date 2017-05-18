#include "config/Base.hh"

#include "tests/Base.hh"

#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include "topology/io/Function.hh"

#include <set>

template <class D, class V> void test( const std::string& filename )
{
  ALEPH_TEST_BEGIN( "Functions file parsing" );

  using Simplex           = aleph::topology::Simplex<D, V>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  auto complexes
    = aleph::topology::io::loadFunctions<SimplicialComplex>( filename );

  ALEPH_ASSERT_EQUAL( complexes.size(), 2 );

  auto K = complexes.at(0);
  auto L = complexes.at(1);

  ALEPH_ASSERT_EQUAL( K.size(), L.size() );
  ALEPH_ASSERT_THROW( K == L ); // modulo weights, both complexes should contain
                                // the same simplices

  ALEPH_TEST_END();
}

int main()
{
  std::vector<std::string> inputs = {
    CMAKE_SOURCE_DIR + std::string( "/tests/input/Functions_simple.txt" ),
    CMAKE_SOURCE_DIR + std::string( "/tests/input/Functions_Reeb.txt" )
  };

  for( auto&& input : inputs )
  {
    test<double,unsigned>      ( input );
    test<double,unsigned short>( input );
    test<float, unsigned>      ( input );
    test<float, unsigned short>( input );
  }
}
