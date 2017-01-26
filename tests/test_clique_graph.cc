#include "tests/Base.hh"

#include "topology/CliqueGraph.hh"
#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include <vector>

using namespace aleph::topology;
using namespace aleph;

template <class Data, class Vertex> void triangle()
{
  ALEPH_TEST_BEGIN( "Triangle" );

  using Simplex           = Simplex<Data, Vertex>;
  using SimplicialComplex = SimplicialComplex<Simplex>;

  std::vector<Simplex> simplices
    = { {0}, {1}, {2}, {0,1}, {0,2}, {1,2}, {0,1,2} };

  SimplicialComplex K( simplices.begin(), simplices.end() );

  auto C = getCliqueGraph( K, 1 );

  ALEPH_ASSERT_THROW( C.empty() == false );

  std::cerr << C << "\n";

  ALEPH_TEST_END();
}

int main()
{
  triangle<double, unsigned>();
}
