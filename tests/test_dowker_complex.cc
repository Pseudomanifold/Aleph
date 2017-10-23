#include <tests/Base.hh>

#include <aleph/geometry/DowkerComplex.hh>

#include <vector>

template <class T> void test()
{
  ALEPH_TEST_BEGIN( "Simple directed networks" );

  using Matrix = std::vector< std::vector<T> >;

  Matrix X( 3, std::vector<T>( 3 ) );

  X[0][1] = T(6);
  X[0][2] = T(4);
  X[1][0] = T(1);
  X[1][2] = T(5);
  X[2][0] = T(2);
  X[2][1] = T(3);

  Matrix Y( 3, std::vector<T>( 3 ) );

  Y[0][1] = T(6);
  Y[0][2] = T(2);
  Y[1][0] = T(1);
  Y[1][2] = T(5);
  Y[2][0] = T(4);
  Y[2][1] = T(3);

  auto X_pairs = aleph::geometry::admissiblePairs( X, T(6) );
  auto Y_pairs = aleph::geometry::admissiblePairs( Y, T(6) );

  ALEPH_ASSERT_THROW( X_pairs.empty() == false );
  ALEPH_ASSERT_THROW( Y_pairs.empty() == false );

  ALEPH_ASSERT_EQUAL( X.size(), Y.size() );

  ALEPH_TEST_END();
}

int main( int, char** )
{
  test<float> ();
  test<double>();
}
