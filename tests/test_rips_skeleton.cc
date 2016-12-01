#include "complexes/FLANN.hh"
#include "complexes/RipsSkeleton.hh"

#include "config/Base.hh"

#include "containers/PointCloud.hh"

#include "distances/Euclidean.hh"

#include "tests/Base.hh"

#include <algorithm>
#include <vector>

using namespace aleph::complexes;
using namespace aleph;

template <class T> void test()
{
  ALEPH_TEST_BEGIN( "Rips skeleton test with different types" );

  using PointCloud = PointCloud<T>;
  using Distance   = aleph::distances::Euclidean<T>;

  PointCloud pointCloud = load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_colon_separated.txt" ) );

  ALEPH_ASSERT_THROW( pointCloud.size()      == 150 );
  ALEPH_ASSERT_THROW( pointCloud.dimension() ==   4);

  using FLANN = FLANN<PointCloud, Distance>;

  FLANN flannWrapper( pointCloud );
  RipsSkeleton<FLANN> ripsSkeleton;

  auto K        = ripsSkeleton.build( flannWrapper, 8.0 );
  auto numEdges = std::count_if( K.begin(), K.end(),
                                     [] ( const typename decltype(K)::ValueType& s )
                                     {
                                       return s.dimension() == 1;
                                     } );

  ALEPH_ASSERT_THROW( K.empty() == false );
  ALEPH_ASSERT_THROW( numEdges > 0 );
  ALEPH_ASSERT_THROW( numEdges == ( pointCloud.size() * ( pointCloud.size() - 1 ) ) / 2 );

  ALEPH_TEST_END();
}

int main()
{
  test<float> ();
  test<double>();
}
