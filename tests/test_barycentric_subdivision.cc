#include <tests/Base.hh>

#include <aleph/topology/BarycentricSubdivision.hh>
#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

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

template <class T> void testWeighted()
{
  ALEPH_TEST_BEGIN( "Weighted triangle" );

  using DataType          = float;
  using VertexType        = T;
  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K = {
   {0},
   {1},
   {2},
   Simplex( {0,1},   DataType(1) ), Simplex( {0,2}, DataType(2) ),
   Simplex( {1,2},   DataType(1) ),
   Simplex( {0,1,2}, DataType(2) )
  };

  aleph::topology::BarycentricSubdivision Sd;

  auto L = Sd(K, [] ( std::size_t dimension ) { return dimension == 0 ? 0 : 0.5; } );

  {
    bool useMaximum                  = true;
    bool skipOneDimensionalSimplices = true;

    L.recalculateWeights( useMaximum, skipOneDimensionalSimplices );
    L.sort( aleph::topology::filtrations::Data<Simplex>() );
  }

  ALEPH_ASSERT_THROW( L.empty() == false );
  ALEPH_ASSERT_THROW( K.size() < L.size() );
  ALEPH_ASSERT_EQUAL( L.size(), 25 );

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

  testWeighted<short>   ();
  testWeighted<int>     ();
  testWeighted<unsigned>();
  testWeighted<long>    ();
}
