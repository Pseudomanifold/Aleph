#include <aleph/config/Base.hh>
#include <aleph/config/FLANN.hh>

#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/FLANN.hh>
#include <aleph/geometry/RipsSkeleton.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <tests/Base.hh>

#include <algorithm>
#include <vector>

using namespace aleph::containers;
using namespace aleph::geometry;
using namespace aleph;

template <class T> void test()
{
  ALEPH_TEST_BEGIN( "Rips skeleton test with different types" );

  using PointCloud = PointCloud<T>;
  using Distance   = distances::Euclidean<T>;

  PointCloud pointCloud = load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_colon_separated.txt" ) );

  ALEPH_ASSERT_THROW( pointCloud.size()      == 150 );
  ALEPH_ASSERT_THROW( pointCloud.dimension() ==   4);

#ifdef ALEPH_WITH_FLANN
  using Wrapper = FLANN<PointCloud, Distance>;
#else
  using Wrapper = BruteForce<PointCloud, Distance>;
#endif

  Wrapper wrapper( pointCloud );
  RipsSkeleton<Wrapper> ripsSkeleton;

  auto K        = ripsSkeleton( wrapper, 8.0 );
  auto numEdges = std::count_if( K.begin(), K.end(),
                                     [] ( const typename decltype(K)::ValueType& s )
                                     {
                                       return s.dimension() == 1;
                                     } );

  ALEPH_ASSERT_THROW( K.empty() == false );
  ALEPH_ASSERT_THROW( numEdges > 0 );
  ALEPH_ASSERT_THROW( static_cast<std::size_t>( numEdges ) == ( pointCloud.size() * ( pointCloud.size() - 1 ) ) / 2 );

  ALEPH_TEST_END();
}

int main()
{
  test<float> ();
  test<double>();
}
