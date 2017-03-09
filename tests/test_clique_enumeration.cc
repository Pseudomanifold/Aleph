#include "tests/Base.hh"

#include "topology/MaximalCliques.hh"

#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include <algorithm>
#include <vector>

using namespace aleph::topology;
using namespace aleph;

template <class Data, class Vertex> void triangles()
{
  ALEPH_TEST_BEGIN( "Triangles [connected & unconnected]" );

  using Simplex           = Simplex<Data, Vertex>;
  using SimplicialComplex = SimplicialComplex<Simplex>;

  // 2---1
  // |  /|
  // | / |
  // |/  |
  // 0---3
  //
  // Expected cliques: {0,1,2}, {0,3,1}
  std::vector<Simplex> trianglesConnected
    = {
        {0,1}, {0,2}, {0,3}, {1,2}, {1,3},
        {0,1,2}, {0,1,3}
    };

  // 2---1   5
  // |  /   /|
  // | /   / |
  // |/   /  |
  // 0---3---4
  //
  // Expected cliques: {0,3}, {0,1,2}, {3,4,5}
  std::vector<Simplex> trianglesDisconnected
    = {
        {0,1}, {0,2}, {0,3}, {1,2}, {3,4}, {3,5}, {4,5},
        {0,1,2}, {3,4,5}
    };

  SimplicialComplex K1( trianglesConnected.begin()   , trianglesConnected.end() );
  SimplicialComplex K2( trianglesDisconnected.begin(), trianglesDisconnected.end() );

  auto C11 = maximalCliquesBronKerbosch( K1 );
  auto C12 = maximalCliquesKoch( K1 );
  auto C21 = maximalCliquesBronKerbosch( K2 );
  auto C22 = maximalCliquesKoch( K2 );

  ALEPH_ASSERT_THROW( C11.empty() == false );
  ALEPH_ASSERT_THROW( C12.empty() == false );
  ALEPH_ASSERT_THROW( C21.empty() == false );
  ALEPH_ASSERT_THROW( C22.empty() == false );

  ALEPH_ASSERT_EQUAL( C11.size(), C12.size() );
  ALEPH_ASSERT_EQUAL( C21.size(), C22.size() );

  ALEPH_TEST_END();
}

int main()
{
  triangles<double, unsigned>();
  triangles<float,  unsigned>();
}
