#include <tests/Base.hh>

#include <aleph/math/KernelDensityEstimator.hh>

#include <vector>

void test1D()
{
  ALEPH_TEST_BEGIN( "1D example" );

  std::vector<double> data      = { -2.1, -1.3, -0.4, 1.9, 5.1, 6.2 };
  std::vector<double> densities = { 0.1074, 0.1244, 0.1184, 0.0691, 0.0828, 0.0789 };

  aleph::math::KernelDensityEstimator kde( 1, 1 );

  for( std::size_t i = 0; i < data.size(); i++ )
  {
    double density = kde( data.begin(), data.end(),
                          data[i],
                          aleph::math::kernels::Gaussian( std::sqrt( 2.25 ) ),
                          aleph::math::norms::Identity() );

    ALEPH_ASSERT_THROW( std::abs( densities[i] - density ) < 1e-4 );
  }

  ALEPH_TEST_END();
}

void test2D()
{
}

int main(int, char**)
{
  test1D();
  test2D();
}
