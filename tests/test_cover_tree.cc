#include <tests/Base.hh>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

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
  ALEPH_ASSERT_THROW( ct.checkCoveringInvariant() );
  ALEPH_ASSERT_THROW( ct.checkSeparatingInvariant() );

  ALEPH_TEST_END();
}

template <class T> void testSimplePermutations()
{
  ALEPH_TEST_BEGIN( "Simple (using permutations)" );

  std::vector<T> data = {7,8,9,10,11,12,13};

  do
  {
    CoverTree<T,
              SimpleMetric<T> > ct;

    // Debug output ----------------------------------------------------

    std::cerr << "Permutation: ";

    for( auto&& x : data )
      std::cerr << x << " ";

    std::cerr << "\n";

    // Check validity of tree ------------------------------------------

    for( auto&& x : data )
      ct.insert( x );

    ALEPH_ASSERT_THROW( ct.checkLevelInvariant() );
    ALEPH_ASSERT_THROW( ct.checkCoveringInvariant() );
    ALEPH_ASSERT_THROW( ct.checkSeparatingInvariant() );
  }
  while( std::next_permutation( data.begin(), data.end() ) );

  ALEPH_TEST_END();
}

template <class T> struct Point
{
  T x;
  T y;

  bool operator<( const Point& other ) const noexcept
  {
    if( x == other.x )
      return y < other.y;
    else
      return x < other.x;
  }
};

template <class T> std::ostream& operator<<( std::ostream& o, const Point<T>& p )
{
  o << p.x << "," << p.y;
  return o;
}

template <class T> struct EuclideanMetric
{
  T operator()( Point<T> a, Point<T> b )
  {
    return std::sqrt( std::pow( a.x - b.x, T(2) ) + std::pow( a.y - b.y, T(2) ) );
  }
};

template <class T> bool contains( const Point<T>& centre, const Point<T>& p, T r )
{
  EuclideanMetric<T> metric;
  return metric( centre, p ) <= r;
}

template <class T> void test2D()
{
  ALEPH_TEST_BEGIN( "2D" );

  using Point     = Point<T>;
  using Metric    = EuclideanMetric<T>;
  using CoverTree = CoverTree<Point, Metric>;

  CoverTree ct;

  std::ifstream in( CMAKE_SOURCE_DIR + std::string( "/tests/input/Cover_tree_simple.txt" ) );
  ALEPH_ASSERT_THROW( in );

  std::vector<Point> points;

  std::string line;
  while( std::getline( in, line ) )
  {
    std::stringstream converter( line );

    T x = T();
    T y = T();

    converter >> x
              >> y;

    ALEPH_ASSERT_THROW( not converter.fail() );

    points.push_back( {x,y} );
  }

  ALEPH_ASSERT_EQUAL( points.size(), 15 );

  for( auto&& p : points )
    ct.insert( p );

  ALEPH_ASSERT_THROW( ct.isValid() );

  auto&& nodesByLevel = ct.getNodesByLevel();

  // Determine radii, i.e. *level* of the original data set. Afterwards,
  // using the corresponding point as the centre, we can check how often
  // certain points are being covered.

  std::map<Point, unsigned> covered;

  for( auto&& pair : nodesByLevel )
  {
    auto&& level  = pair.first;
    auto&& centre = pair.second;

    for( auto&& p : points )
    {
      // TODO: fix radius/level calculation; is this an implementation
      // detail of the tree?
      if( contains( centre, p, T( std::pow( T(2), level ) ) ) )
        covered[p] += 1;
    }
  }

  for( auto&& pair : covered )
    std::cerr << pair.first << ": " << pair.second << "\n";

  ALEPH_ASSERT_EQUAL( nodesByLevel.size(), points.size() );
  ALEPH_TEST_END();
}

int main( int, char** )
{
  testSimple<double>();
  testSimple<float> ();

  testSimplePermutations<double>();
  testSimplePermutations<float> ();

  test2D<double>();
  test2D<float> ();
}
