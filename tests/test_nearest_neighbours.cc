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

  using IndexType   = typename FLANN<PointCloud>::IndexType;
  using ElementType = typename FLANN<PointCloud>::ElementType;

  std::vector< std::vector<IndexType> > indices;
  std::vector< std::vector<ElementType> > distances;

  flannWrapper.radiusSearch( static_cast<T>( 0.5 ), indices, distances );

  for( auto&& i : indices )
    for( auto&& j : i )
      std::cerr << j << " ";

  ALEPH_TEST_END();
}

int main()
{
  test<float> ();
  test<double>();
}
