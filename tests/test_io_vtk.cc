#include "config/Base.hh"

#include "tests/Base.hh"

#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include "topology/io/VTK.hh"

template <class D, class V> void test()
{
  ALEPH_TEST_BEGIN( "VTK structured grid parsing" );

  using Simplex           = aleph::topology::Simplex<D, V>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K;

  aleph::topology::io::VTKStructuredGridReader reader;
  reader( CMAKE_SOURCE_DIR + std::string( "/tests/input/Simple.vtk" ),
          K );

  ALEPH_TEST_END();
}

int main()
{
  test<double,unsigned>      ();
  test<double,unsigned short>();
  test<float, unsigned>      ();
  test<float, unsigned short>();
}
