#include <tests/Base.hh>

#include <aleph/topology/BarycentricSubdivision.hh>
#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <vector>

template <class T> void test()
{
  using DataType          = float;
  using VertexType        = T;
  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K = {
   {0},
   {1},
   {2},
   {0,1}, {0,2},
   {1,2},
   {0,1,2}
  };

  aleph::topology::BarycentricSubdivision Sd;

  auto L = Sd(K);

  ALEPH_ASSERT_THROW( L.empty() == false );
  ALEPH_ASSERT_THROW( K.size() < L.size() );
}

int main(int, char**)
{
  test<short>   ();
  test<int>     ();
  test<unsigned>();
  test<long>    ();
}
