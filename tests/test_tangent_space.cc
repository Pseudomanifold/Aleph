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
  auto n = 100;
  auto k = 5;

  // FIXME: use for tangent space estimation
  (void) k;

  using PointCloud = PointCloud<T>;
  PointCloud pc( n, 2 );

  for( int i = 0; i < n; i++ )
  {
    auto phi = T(2*M_PI*i/n);
    auto r   = T(1.0);
    auto x   = r * std::cos( phi );
    auto y   = r * std::sin( phi );

    pc.set(i, {x,y} );
  }

#ifdef ALEPH_WITH_EIGEN
  TangentSpace ts;
  ts( pc );
#endif
}

template <class T> void test()
{
#ifdef ALEPH_WITH_EIGEN
  using PointCloud = PointCloud<T>;

  PointCloud pc = load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_colon_separated.txt" ) );

  TangentSpace ts;
  ts( pc );
#endif
}

int main( int, char** )
{
  testCircle<float> ();
  testCircle<double>();

  //test<float> ();
  //test<double>();
}
