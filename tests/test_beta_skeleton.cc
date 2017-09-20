#include <aleph/config/Base.hh>

#include <tests/Base.hh>

#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/BetaSkeleton.hh>

#include <aleph/geometry/distances/Euclidean.hh>

using namespace aleph;
using namespace containers;
using namespace geometry;
using namespace distances;

template <class T> void test()
{
  using PointCloud = PointCloud<T>;
  using Distance   = Euclidean<T>;

  PointCloud pc = load<T>( CMAKE_SOURCE_DIR + std::string("/tests/input/Iris_tab_separated.txt") );

  auto K0 = buildBetaSkeletonNaive( pc, 0, Distance() );
  auto K1 = buildBetaSkeletonNaive( pc, 1, Distance() );
  auto K2 = buildBetaSkeletonNaive( pc, 2, Distance() );
  auto K3 = buildBetaSkeletonNaive( pc, 3, Distance() );

  for( auto&& K : {K0,K1,K2,K3} )
    ALEPH_ASSERT_THROW( K.empty() == false );
}

int main( int, char** )
{
  test<float> ();
  test<double>();
}
