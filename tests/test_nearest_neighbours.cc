#include "complexes/FLANN.hh"
#include "complexes/NearestNeighbours.hh"

#include "config/Base.hh"

#include "containers/PointCloud.hh"

#include "distances/Euclidean.hh"

#include "tests/Base.hh"

#include <vector>

#include <cassert>

using namespace aleph::complexes;
using namespace aleph;

template <class T> void test()
{
  ALEPH_TEST_BEGIN( "Nearest-neighbour calculation with different types" );

  using PointCloud = PointCloud<T>;
  using Distance   = aleph::distances::Euclidean<T>;

  PointCloud pointCloud = load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_colon_separated.txt" ) );

  ALEPH_ASSERT_THROW( pointCloud.size()      == 150 );
  ALEPH_ASSERT_THROW( pointCloud.dimension() ==   4);

  FLANN<PointCloud, Distance> flannWrapper( pointCloud );

  using IndexType   = typename FLANN<PointCloud, Distance>::IndexType;
  using ElementType = typename FLANN<PointCloud, Distance>::ElementType;

  std::vector< std::vector<IndexType> > indices;
  std::vector< std::vector<ElementType> > distances;

  // Check that an *empty* radius does not return any indices

  flannWrapper.radiusSearch( static_cast<T>( 0.0 ), indices, distances );

  ALEPH_ASSERT_THROW( indices.size() == pointCloud.size() );
  for( auto&& i : indices )
    ALEPH_ASSERT_THROW( i.empty() == true );

  // Check that a large radius returns *all* indices

  flannWrapper.radiusSearch( static_cast<T>( 8.0 ), indices, distances );

  ALEPH_ASSERT_THROW( indices.size() == pointCloud.size() );
  for( auto&& i : indices )
    ALEPH_ASSERT_THROW( i.size() == pointCloud.size() );

  ALEPH_TEST_END();
}

int main()
{
  test<float> ();
  test<double>();
}
