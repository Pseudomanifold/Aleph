#include <aleph/containers/FractalDimension.hh>
#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <tests/Base.hh>

using namespace aleph::geometry::distances;
using namespace aleph::containers;
using namespace aleph;

template <class T> void testCorrelationDimension()
{
  ALEPH_TEST_BEGIN( "Correlation dimension" );

  using PointCloud = PointCloud<T>;
  PointCloud pc = load<T>(
    CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_comma_separated.txt" )
  );

  auto cds = correlationDimensionIntegral( pc, Euclidean<double>() );

  ALEPH_ASSERT_THROW( cds.x.empty() == false );
  ALEPH_ASSERT_THROW( cds.y.empty() == false );

  auto nu = correlationDimension( cds );

  ALEPH_ASSERT_THROW( nu > 0.0 );
  ALEPH_ASSERT_THROW( nu > 1.0 );

  ALEPH_TEST_END();
}

int main(int, char**)
{
  testCorrelationDimension<float> ();
  testCorrelationDimension<double>();
}
