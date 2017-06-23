#include <aleph/config/Base.hh>

#include <aleph/persistentHomology/ConnectedComponents.hh>

#include <tests/Base.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/topology/io/VTK.hh>

#include <algorithm>

template <class D, class V> void test()
{
  ALEPH_TEST_BEGIN( "VTK structured grid parsing" );

  using Simplex           = aleph::topology::Simplex<D, V>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K;

  aleph::topology::io::VTKStructuredGridReader reader;
  reader( CMAKE_SOURCE_DIR + std::string( "/tests/input/Simple.vtk" ),
          K,
          [] ( D a, D b ) { return std::min(a,b); } );

  K.sort( aleph::topology::filtrations::Data<Simplex, std::greater<D> >() );

  auto n0 = std::count_if( K.begin(), K.end(), [] ( const Simplex& s ) { return s.dimension() == 0; } );
  auto n1 = std::count_if( K.begin(), K.end(), [] ( const Simplex& s ) { return s.dimension() == 1; } );
  auto n2 = std::count_if( K.begin(), K.end(), [] ( const Simplex& s ) { return s.dimension() == 2; } );

  ALEPH_ASSERT_EQUAL( n0,  5000 );
  ALEPH_ASSERT_EQUAL( n1, 12300 ); // TODO: check whether this is correct
  ALEPH_ASSERT_EQUAL( n2,     0 );

  // TODO:
  //  - Count number of 'regular' vertices
  //  - Count number of 'irregular' vertices (392 = 2*(2*nx+2*ny-4))

  auto pdpp = aleph::calculateZeroDimensionalPersistenceDiagram( K );
  auto pd   = std::get<0>( pdpp );

  ALEPH_ASSERT_EQUAL( pd.size(), 3 );

  ALEPH_TEST_END();
}

int main()
{
  test<double,unsigned>      ();
  test<double,unsigned short>();
  test<float, unsigned>      ();
  test<float, unsigned short>();
}
