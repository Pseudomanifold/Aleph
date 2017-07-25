#include <tests/Base.hh>

#include <aleph/topology/BarycentricSubdivision.hh>
#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <algorithm>
#include <vector>

template <class T> void test()
{
  ALEPH_TEST_BEGIN( "Triangle" );

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

  auto num0Simplices = std::count_if( L.begin(), L.end(), [] ( const Simplex& s ) { return s.dimension() == 0; } );
  auto num1Simplices = std::count_if( L.begin(), L.end(), [] ( const Simplex& s ) { return s.dimension() == 1; } );
  auto num2Simplices = std::count_if( L.begin(), L.end(), [] ( const Simplex& s ) { return s.dimension() == 2; } );

  ALEPH_ASSERT_EQUAL( num0Simplices,  7 );
  ALEPH_ASSERT_EQUAL( num1Simplices, 12 );
  ALEPH_ASSERT_EQUAL( num2Simplices,  6 );

  ALEPH_TEST_END();
}

int main(int, char**)
{
  test<short>   ();
  test<int>     ();
  test<unsigned>();
  test<long>    ();
}
