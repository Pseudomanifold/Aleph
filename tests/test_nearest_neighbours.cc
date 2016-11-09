#include "complexes/FLANN.hh"
#include "complexes/NearestNeighbours.hh"

#include "config/Base.hh"

#include "containers/PointCloud.hh"

#include "tests/Base.hh"

#include <vector>

#include <cassert>

using namespace aleph::complexes;
using namespace aleph;

template <class T> void test()
{
  ALEPH_TEST_BEGIN( "Nearest-neighbour calculation with different types" );

  using PointCloud = PointCloud<T>;

  PointCloud pointCloud = load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_colon_separated.txt" ) );

  ALEPH_ASSERT_THROW( pointCloud.size()      == 150 );
  ALEPH_ASSERT_THROW( pointCloud.dimension() ==   4);

  FLANN<PointCloud> flannWrapper( pointCloud );

  ALEPH_TEST_END();
}

int main()
{
  test<float> ();
  test<double>();
}
