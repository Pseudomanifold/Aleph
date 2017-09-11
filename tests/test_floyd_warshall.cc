#include <tests/Base.hh>

#include <aleph/topology/FloydWarshall.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <vector>

template <class T> void test()
{
  using Simplex           = aleph::topology::Simplex<T, unsigned>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  std::vector<Simplex> simplices
    = {
    {1}, {2}, {3}, {4},
    Simplex( {1,2}, T(1) ),
    Simplex( {2,3}, T(2) ),
    Simplex( {3,4}, T(3) ),
    Simplex( {4,1}, T(4) ),
    Simplex( {4,2}, T(7) )
  };

  SimplicialComplex K( simplices.begin(), simplices.end() );
  auto M = aleph::topology::floydWarshall( K );

  ALEPH_ASSERT_THROW( M.empty() == false );
  ALEPH_ASSERT_EQUAL( M.numRows(), 4 );
  ALEPH_ASSERT_EQUAL( M(0,2), T(3) );
  ALEPH_ASSERT_EQUAL( M(3,1), T(5) );
  ALEPH_ASSERT_EQUAL( M(0,0), T(0) );
}

int main( int, char** )
{
  test<float> ();
  test<double>();
}
