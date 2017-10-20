#include <tests/Base.hh>

#include <aleph/topology/Partitions.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <vector>

// FIXME: remove after debugging
#include <iostream>

template <class T> void test()
{
  ALEPH_TEST_BEGIN( "Bisection" );

  using Simplex           = aleph::topology::Simplex<T>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;
  using VertexType        = typename Simplex::VertexType;

  SimplicialComplex K;

  K.push_back( Simplex( VertexType(0), 0 ) ); K.push_back( Simplex( VertexType(1), 0 ) );
  K.push_back( Simplex( VertexType(2), 0 ) ); K.push_back( Simplex( VertexType(3), 0 ) );
  K.push_back( Simplex( VertexType(4), 0 ) ); K.push_back( Simplex( VertexType(5), 0 ) );

  K.push_back( Simplex( {0,1}, 1 ) );
  K.push_back( Simplex( {0,2}, 1 ) );
  K.push_back( Simplex( {1,2}, 1 ) );

  K.push_back( Simplex( {3,4}, 1 ) );
  K.push_back( Simplex( {3,5}, 1 ) );
  K.push_back( Simplex( {4,5}, 1 ) );

  K.push_back( Simplex( {2,3}, 1 ) );

  K.push_back( Simplex( {3,4,5}, 1 ) );
  K.push_back( Simplex( {0,1,2}, 1 ) );

  // Technically, this invalidates the simplicial complex, but I am only
  // interested in figuring out whether the partition will *ignore* this
  // simplex.
  K.push_back( Simplex( {0,1,2,3}, 1 ) );

  auto complexes = bisect( K );

  ALEPH_ASSERT_EQUAL( complexes.size(), 2 );

  auto&& K1 = complexes.front();
  auto&& K2 = complexes.back();

  ALEPH_ASSERT_EQUAL( K1.size(), K2.size() );

  ALEPH_TEST_END();
}

int main( int, char** )
{
  test<float> ();
  test<double>();
}
