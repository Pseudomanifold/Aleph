#include <tests/Base.hh>

#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/CechComplex.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <algorithm>
#include <iterator>
#include <set>
#include <vector>

#include <cmath>

using namespace aleph::containers;
using namespace aleph::geometry;
using namespace aleph::topology;
using namespace aleph;

template <class T> void triangle()
{
  ALEPH_TEST_BEGIN( "Triangle" );

  using PointCloud = PointCloud<T>;
  PointCloud pc    = load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Triangle_point_cloud.txt" ) );

  ALEPH_ASSERT_EQUAL( pc.dimension(), 2);
  ALEPH_ASSERT_EQUAL( pc.size()     , 3);

  ALEPH_TEST_END();
}

int main()
{
  triangle<double>();
  triangle<float> ();
}
