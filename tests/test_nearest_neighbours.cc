
#include "config/Base.hh"
#include "config/FLANN.hh"

#include "containers/PointCloud.hh"

#include "distances/Euclidean.hh"

#include "geometry/BruteForce.hh"
#include "geometry/FLANN.hh"
#include "geometry/NearestNeighbours.hh"

#include "tests/Base.hh"

#include <vector>

#include <cassert>

using namespace aleph::geometry;
using namespace aleph;

template <class Wrapper, class PointCloud> void testInternal( const PointCloud& pointCloud )
{
  Wrapper wrapper( pointCloud );

  using IndexType   = typename Wrapper::IndexType;
  using ElementType = typename Wrapper::ElementType;

  std::vector< std::vector<IndexType> > indices;
  std::vector< std::vector<ElementType> > distances;

  // Check that an *empty* radius does not return any indices

  wrapper.radiusSearch( static_cast<ElementType>( 0.0 ), indices, distances );

  ALEPH_ASSERT_THROW( indices.size() == pointCloud.size() );
  for( auto&& i : indices )
    ALEPH_ASSERT_THROW( i.empty() == true );

  // Check that a large radius returns *all* indices

  wrapper.radiusSearch( static_cast<ElementType>( 8.0 ), indices, distances );

  ALEPH_ASSERT_THROW( indices.size() == pointCloud.size() );
  for( auto&& i : indices )
    ALEPH_ASSERT_THROW( i.size() == pointCloud.size() );
}

template <class T> void test()
{
  ALEPH_TEST_BEGIN( "Nearest-neighbour calculation with different types" );

  using PointCloud = PointCloud<T>;
  using Distance   = aleph::distances::Euclidean<T>;

  PointCloud pointCloud = load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_colon_separated.txt" ) );

  ALEPH_ASSERT_THROW( pointCloud.size()      == 150 );
  ALEPH_ASSERT_THROW( pointCloud.dimension() ==   4);

#ifdef ALEPH_WITH_FLANN
  testInternal< FLANN<PointCloud, Distance> >( pointCloud );
#endif
  testInternal< BruteForce<PointCloud, Distance> >( pointCloud );

  ALEPH_TEST_END();
}

int main()
{
  test<float> ();
  test<double>();
}
