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
  // Expected cliques graph: {0,3}, {0,1,2}, {3,4,5}
  std::vector<Simplex> trianglesDisconnected
    = {
        {0,1}, {0,2}, {0,3}, {1,2}, {3,4}, {3,5}, {4,5},
        {0,1,2}, {3,4,5}
    };

  SimplicialComplex K1( trianglesConnected.begin()   , trianglesConnected.end() );
  SimplicialComplex K2( trianglesDisconnected.begin(), trianglesDisconnected.end() );

  auto C1 = maximalCliquesKoch( K1 );
  auto C2 = maximalCliquesKoch( K2 );

  ALEPH_ASSERT_THROW( C1.empty() == false );
  ALEPH_ASSERT_THROW( C2.empty() == false );

  ALEPH_TEST_END();
}

int main()
{
  triangle<double, unsigned>();
  triangle<float,  unsigned>();

  triangles<double, unsigned>();
  triangles<float,  unsigned>();
}
