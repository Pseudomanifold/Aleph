#include <tests/Base.hh>

#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/TangentSpace.hh>

#include <iostream>

#include <cmath>

using namespace aleph;
using namespace containers;
using namespace geometry;

template <class T> void testCircle()
{
  unsigned n = 200;

  using PointCloud = PointCloud<T>;
  PointCloud pc( n, 2 );

  for( unsigned i = 0; i < n; i++ )
  {
    auto phi = T(2*M_PI*i/(n-1));
    auto r   = T(1.0);
    auto x   = r * std::cos( phi );
    auto y   = r * std::sin( phi );

    pc.set(i, {x,y} );
  }

#ifdef ALEPH_WITH_EIGEN
#if EIGEN_VERSION_AT_LEAST(3,3,0)
  TangentSpace ts;

  unsigned k     = 10;
  auto curvature = ts( pc, k );

  ALEPH_ASSERT_EQUAL( curvature.size(), std::size_t(n) );

  // TODO: more tests for curvature here...

  for( auto&& c : curvature )
    std::cout << c << "\n";
#endif
#endif
}

template <class T> void test()
{
#ifdef ALEPH_WITH_EIGEN
#if EIGEN_VERSION_AT_LEAST(3,3,0)
  using PointCloud = PointCloud<T>;

  PointCloud pc = load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_colon_separated.txt" ) );

  TangentSpace ts;
  ts( pc, 10 );
#endif
#endif
}

int main( int, char** )
{
  testCircle<float> ();
  testCircle<double>();

  test<float> ();
  test<double>();
}
