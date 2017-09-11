#include <tests/Base.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/CombinatorialCurvature.hh>

#include <iterator>
#include <vector>

template <class T > void testSphere()
{
  ALEPH_TEST_BEGIN( "Sphere" );

  using DataType   = bool;
  using VertexType = T;

  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K( {
      {0,1,2,3},
      {0,1,2}, {0,1,3}, {1,2,3}, {0,2,3},
      {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3},
      {0}, {1}, {2}, {3}
    }
  );

  K.sort();

  std::vector<T> curvature;
  aleph::topology::curvature( K, std::back_inserter( curvature ) );

  ALEPH_ASSERT_THROW( curvature.empty() == false );
  ALEPH_ASSERT_EQUAL( curvature.size(), 6 );
  ALEPH_ASSERT_EQUAL( curvature.front(), curvature.back() );
  ALEPH_ASSERT_EQUAL( curvature.front(), 4 );

  ALEPH_TEST_END();
}

int main( int, char** )
{
  testSphere<unsigned>      ();
  testSphere<unsigned short>();
}
