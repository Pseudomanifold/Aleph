#include <tests/Base.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/io/HDF5.hh>

template <class D, class V> void test()
{
  ALEPH_TEST_BEGIN( "HDF5 file simple data set parsing" );

  using Simplex           = aleph::topology::Simplex<D, V>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K;
  SimplicialComplex L;

  aleph::topology::io::HDF5SimpleDataSpaceReader reader;
  reader.setDataSetName( "Simple" );

  reader( CMAKE_SOURCE_DIR + std::string( "/tests/input/Simple.h5" ),
          K,
          [] ( D a, D b ) { return std::min(a,b); } );

  reader( CMAKE_SOURCE_DIR + std::string( "/tests/input/Simple.h5" ), L );

  ALEPH_ASSERT_THROW( K.empty() == false );
  ALEPH_ASSERT_THROW( L.empty() == false );

  ALEPH_ASSERT_EQUAL( K.size(), L.size() );

  auto checkSimplexCount = [] ( const SimplicialComplex& K )
  {
    auto n0 = std::count_if( K.begin(), K.end(), [] ( const Simplex& s ) { return s.dimension() == 0; } );
    auto n1 = std::count_if( K.begin(), K.end(), [] ( const Simplex& s ) { return s.dimension() == 1; } );
    auto n2 = std::count_if( K.begin(), K.end(), [] ( const Simplex& s ) { return s.dimension() == 2; } );

    ALEPH_ASSERT_EQUAL( n0,  9 );
    ALEPH_ASSERT_EQUAL( n1, 16 );
    ALEPH_ASSERT_EQUAL( n2,  8 );
  };

  checkSimplexCount( K );
  checkSimplexCount( L );

  ALEPH_TEST_END();
}

int main(int, char**)
{
  test<double,unsigned>      ();
  test<double,unsigned short>();
  test<float, unsigned>      ();
  test<float, unsigned short>();
}
