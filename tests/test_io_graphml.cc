#include <tests/Base.hh>

#include <aleph/config/TinyXML2.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/io/GML.hh>
#include <aleph/topology/io/SimplicialComplexReader.hh>

#include <set>

template <class D, class V> void test( const std::string& filename )
{
  ALEPH_TEST_BEGIN( "GraphML file parsing" );

  using Simplex           = aleph::topology::Simplex<D, V>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K;

  aleph::topology::io::GMLReader reader;
  reader( filename, K );

  ALEPH_ASSERT_THROW( K.empty() == false );

  ALEPH_TEST_END();
}

int main()
{
  auto input = CMAKE_SOURCE_DIR + std::string( "/tests/input/Simple.xml" );

#ifdef ALEPH_WITH_TINYXML2
  test<double,unsigned>      ( input );
  test<double,unsigned short>( input );
  test<float, unsigned>      ( input );
  test<float, unsigned short>( input );
#endif
}
