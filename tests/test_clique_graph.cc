#include "tests/Base.hh"

#include "topology/CliqueGraph.hh"
#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include <algorithm>
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
  ALEPH_ASSERT_THROW( C.size()  == 6 );

  ALEPH_TEST_END();
}

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
  // Expected clique graph: {0,1,2} -- {0,3,1}
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
  // Expected clique graph: {0,1,2} -- {0,3,1}
  std::vector<Simplex> trianglesDisconnected
    = {
        {0,1}, {0,2}, {0,3}, {1,2}, {3,4}, {3,5}, {4,5},
        {0,1,2}, {3,4,5}
    };

  SimplicialComplex K1( trianglesConnected.begin()   , trianglesConnected.end() );
  SimplicialComplex K2( trianglesDisconnected.begin(), trianglesDisconnected.end() );

  auto C1 = getCliqueGraph( K1, 2 );
  auto C2 = getCliqueGraph( K2, 2 );

  ALEPH_ASSERT_THROW( C1.empty() == false );
  ALEPH_ASSERT_THROW( C1.size()  == 3     );
  ALEPH_ASSERT_THROW( std::count_if( C1.begin(), C1.end(), [] ( const Simplex& s ) { return s.dimension() == 0; } ) == 2 );
  ALEPH_ASSERT_THROW( std::count_if( C1.begin(), C1.end(), [] ( const Simplex& s ) { return s.dimension() == 1; } ) == 1 );

  ALEPH_ASSERT_THROW( C2.empty() == false );
  ALEPH_ASSERT_THROW( C2.size()  == 2     );
  ALEPH_ASSERT_THROW( std::count_if( C2.begin(), C2.end(), [] ( const Simplex& s ) { return s.dimension() == 0; } ) == 2 );
  ALEPH_ASSERT_THROW( std::count_if( C2.begin(), C2.end(), [] ( const Simplex& s ) { return s.dimension() == 1; } ) == 0 );

  ALEPH_TEST_END();
}

int main()
{
  triangle<double, unsigned>();
  triangle<float,  unsigned>();

  triangles<double, unsigned>();
  triangles<float,  unsigned>();
}
