#include <tests/Base.hh>

#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/distances/Euclidean.hh>
#include <aleph/geometry/WitnessComplex.hh>

#include <algorithm>
#include <iterator>
#include <set>

template <class T> void test()
{
  using Distance   = aleph::distances::Euclidean<T>;
  using PointCloud = aleph::containers::PointCloud<T>;

  PointCloud pc( 8, 2 );

  pc.set(0, {-1.0, 0.0} );
  pc.set(1, { 0.0,-1.0} );
  pc.set(2, { 1.0, 0.0} );
  pc.set(3, { 2.0, 1.0} );
  pc.set(4, { 1.0, 1.0} );
  pc.set(5, { 0.0, 2.0} );
  pc.set(6, {-1.0, 1.0} );
  pc.set(7, {-2.0, 1.0} );

  std::vector<std::size_t> indices = {0,2,4,6};

  auto K
    = aleph::geometry::buildWitnessComplex<Distance>(
        pc, indices.begin(), indices.end() );

  using Simplex    = typename decltype(K)::ValueType;
  using VertexType = typename Simplex::VertexType;

  {
    std::set<VertexType> vertices;

    K.vertices( std::inserter( vertices, vertices.begin() ) );

    ALEPH_ASSERT_EQUAL( vertices.size(), indices.size() );
  }

  auto numEdges = std::count_if( K.begin(), K.end(), [] ( const Simplex& s ) { return s.dimension() == 1; } );

  ALEPH_ASSERT_EQUAL( numEdges, 4 );
}

int main(int, char**)
{
  test<float> ();
  test<double>();
}
