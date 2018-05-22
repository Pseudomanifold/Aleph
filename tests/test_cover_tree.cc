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

  // FIXME
  //std::vector<T> data = {7,8,9,10,11,12,13};

  CoverTree<T,
            SimpleMetric<T> > ct;

  ct.insert( 7 );
  ct.insert( 13 );
  ct.insert( 10 );
  ct.insert( 8 );
  ct.insert( 9 );
  ct.insert( 11 );
  ct.insert( 12 );

  ct.print( std::cerr );

  ALEPH_ASSERT_THROW( ct.checkLevelInvariant() );

  ALEPH_TEST_END();
}

int main( int, char** )
{
  testSimple<double>();
  testSimple<float> ();
}
