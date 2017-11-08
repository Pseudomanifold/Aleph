#include <tests/Base.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>
#include <aleph/topology/Spine.hh>

#include <vector>

template <class T> void testTriangle()
{
  using DataType   = bool;
  using VertexType = T;

  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K = {
    {0,1,2},
    {0,1}, {0,2}, {1,2},
    {0}, {1}, {2}
  };

  auto L = aleph::topology::spine( K );

  ALEPH_ASSERT_THROW( L.size() < K.size() );
  ALEPH_ASSERT_EQUAL( L.size(), 1 );
}

int main( int, char** )
{
  testTriangle<short>   ();
  testTriangle<unsigned>();
}
