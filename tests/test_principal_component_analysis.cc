#include <tests/Base.hh>

#include <aleph/math/PrincipalComponentAnalysis.hh>

#include <algorithm>
#include <iterator>
#include <vector>

#include <cmath>

template <class T> void testSimpleMatrix()
{
  ALEPH_TEST_BEGIN( "Simple matrix" );

  std::vector< std::vector<T> > matrix( 3, std::vector<T>( 2 ) );

  matrix[0][0] = T( 2.0/3.0); matrix[0][1] = T(-3.0 - 2.0/3.0);
  matrix[1][0] = T( 2.0/3.0); matrix[1][1] = T( 4.0 + 1.0/3.0);
  matrix[2][0] = T(-4.0/3.0); matrix[2][1] = T(-2.0/3.0      );

  aleph::math::PrincipalComponentAnalysis pca;
  auto result = pca( matrix );

  std::vector<T> eigenvalues;
  std::transform( result.singularValues.begin(), result.singularValues.end(),
                  std::back_inserter( eigenvalues ),
                  [] ( T value ) { return value * value; } );

  ALEPH_ASSERT_THROW( eigenvalues.empty() == false );
  ALEPH_ASSERT_EQUAL( eigenvalues.size(), 2 );

  ALEPH_ASSERT_THROW( std::abs( eigenvalues[0] - 16.3629 ) < 1e-4 );
  ALEPH_ASSERT_THROW( std::abs( eigenvalues[1] -  1.3037 ) < 1e-4 );

  ALEPH_TEST_END();
}

int main(int, char**)
{
  testSimpleMatrix<float> ();
  testSimpleMatrix<double>();
}
