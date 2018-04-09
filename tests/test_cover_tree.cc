#include <tests/Base.hh>

#include <cmath>

#include <aleph/geometry/CoverTree.hh>

using namespace aleph::geometry;

template <class T> struct SimpleMetric
{
  T operator()( T a, T b )
  {
    return std::abs( a - b );
  }
};

template <class T> void testSimple()
{
  ALEPH_TEST_BEGIN( "Simple" );

  std::vector<T> data = {7,8,9,10,11,12,13};

  CoverTree<T,
            SimpleMetric<T> > ct( data.begin(), data.end() );

  ALEPH_TEST_END();
}

int main( int, char** )
{
  testSimple<double>();
  testSimple<float> ();
}
