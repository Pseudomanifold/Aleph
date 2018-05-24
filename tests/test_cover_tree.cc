#include <tests/Base.hh>

#include <fstream>
#include <iostream>
#include <set>
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

  bool operator==( const Point& other ) const noexcept
  {
    return x == other.x && y == other.y;
  }

  bool operator!=( const Point& other ) const noexcept
  {
    return !this->operator==( other );
  }
};

template <class T> std::ostream& operator<<( std::ostream& o, const Point<T>& p )
{
  o << p.x << " " << p.y;
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

template <class T> T distance( const Point<T>& centre, const Point<T>& p )
{
  EuclideanMetric<T> metric;
  return metric( centre, p );
}

template <class T> std::vector<T> eccentricity( const std::vector< Point<T> >& points )
{
  std::vector<T> E;

  for( auto&& p : points )
  {
    T e = T();

    for( auto&& q : points )
      e += distance( p, q );

    e /= static_cast<T>( points.size() );
    E.push_back( e );
  }

  return E;
}

template <class T> void test2D()
{
  ALEPH_TEST_BEGIN( "2D" );

  using Point     = Point<T>;
  using Metric    = EuclideanMetric<T>;
  using CoverTree = CoverTree<Point, Metric>;

  CoverTree ct;

  std::ifstream in( CMAKE_SOURCE_DIR + std::string( "/tests/input/Cover_tree_sorted.txt" ) );
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

#if 0
  // Check eccentricity-based sorting in order to obtain cover trees
  // with improved balance properties.
  {
    auto e = eccentricity( points );

    std::multimap<T, Point> sortedPoints;

    for( std::size_t i = 0; i < points.size(); i++ )
      sortedPoints.insert( std::make_pair( e.at(i), points.at(i) ) );

    for( auto it = sortedPoints.rbegin(); it != sortedPoints.rend(); ++it )
      ct.insert( it->second );
  }
#endif

  for( auto&& p : points )
    ct.insert( p );

  ALEPH_ASSERT_THROW( ct.isValid() );

  auto&& nodesByLevel = ct.getNodesByLevel();

  // Determine radii, i.e. *level* of the original data set. Afterwards,
  // using the corresponding point as the centre, we can check how often
  // certain points are being covered.

  std::map<Point, unsigned > covered;
  std::multimap<Point, long> levels;
  std::multimap<Point, T   > distances;

  for( auto&& pair : nodesByLevel )
  {
    auto&& level  = pair.first;
    auto&& centre = pair.second;

    for( auto&& p : points )
    {
      // TODO: fix radius/level calculation; is this an implementation
      // detail of the tree?
      if( contains( centre, p, T( std::pow( T(2), level ) ) ) )
      {
        covered[p] += 1;
        levels.insert( std::make_pair( p, level ) );
        distances.insert( std::make_pair( p, distance( centre, p ) ) );
      }
    }
  }

  std::cerr << "# Cover counter\n";

  for( auto&& pair : covered )
    std::cerr << pair.first << ": " << pair.second << "\n";

  std::cerr << "# Levels counter\n";

  for( auto&& p : points )
  {
    std::cerr << p << ": ";

    auto range = levels.equal_range( p );
    for( auto it = range.first; it != range.second; ++it )
     std::cerr << it->second << " ";

    std::cerr << "\n";
  }

  std::cerr << "# Distances counter\n";

  for( auto&& p : points )
  {
    std::cerr << p << ": ";

    auto range = distances.equal_range( p );
    for( auto it = range.first; it != range.second; ++it )
     std::cerr << it->second << " ";

    std::cerr << "\n";
  }

  std::cerr << "# Basic cover distance density\n";

  for( auto&& p : points )
  {
    std::cerr << p << ": ";

    auto range      = distances.equal_range( p );
    double distance = 0.0;

    for( auto it = range.first; it != range.second; ++it )
      distance += it->second;

    distance /= static_cast<T>( std::distance( range.first, range.second ) );

    std::cerr << distance << "\n";
  }

  // DEBUG: output of cover radii --------------------------------------

  {
    std::ofstream out( "/tmp/C.txt" );

    for( auto&& pair : nodesByLevel )
    {
      auto&& level  = pair.first;
      auto&& centre = pair.second;
      auto r        = std::pow( T(2), level );

      out << centre << " " << r << "\n";
    }
  }

  // DEBUG: output of edges --------------------------------------------

  {
    std::set< std::pair<Point, Point> > edges;

    for( auto&& pair : nodesByLevel )
    {
      auto&& level  = pair.first;
      auto&& centre = pair.second;

      for( auto&& p : points )
      {
        // TODO: fix radius/level calculation; is this an implementation
        // detail of the tree?
        if( contains( centre, p, T( std::pow( T(2), level ) ) ) )
        {
          // Induce basic ordering of edges in order to make it easier
          // to print them later on.
          if( centre < p )
            edges.insert( std::make_pair( centre, p ) );
          else
            edges.insert( std::make_pair( p, centre ) );
        }
      }
    }

    std::ofstream out( "/tmp/E.txt" );

    for( auto&& edge : edges )
    {
      out << edge.first  << "\n"
          << edge.second << "\n\n";
    }

    std::cerr << "# Cover tree\n";

    ct.print( std::cerr );

    std::set< std::pair<Point, Point> > filteredEdges;

    auto&& nodesToLevel = ct.nodesToLevel();
    for( auto&& edge : edges )
    {
      auto&& source = edge.first;                 // source node
      auto&& upper  = nodesToLevel.at( source );  // source level
      auto&& target = edge.second;                // target node
      auto&& lower  = nodesToLevel.at( target );  // target level

      auto d       = distance( source, target );
      auto d_lower = std::pow( T(2), lower );
      auto d_upper = std::pow( T(2), upper );

      std::cerr << edge.first << " -- " << edge.second << ":\n"
                << "  " << upper << "," << lower << "," << std::abs( upper - lower ) << "\n";

      if( d <= d_lower && d <= d_upper )
        filteredEdges.insert( edge );
      else
      {
        auto l     = std::min( lower, upper );
        auto L     = std::max( lower, upper );
        unsigned c = 0;

        for( long level = l; level <= L; level++ )
        {
          auto D = std::pow( T(2), level );
          if( d > D )
            c++;
        }

        // TODO: make configurable
        if( c == 1 )
          filteredEdges.insert( edge );
      }

      // TODO: old criterion; this can probably be removed at some point
      // in the future
      #if 0
      else if( d <= d_lower && d > d_upper )
      {
        // Only one criterion holds
      }
      else if( d <= d_upper && d > d_lower )
      {
        // Ditto.
      }
      else
        throw std::runtime_error( "This should never happen" );
      #endif
    }

    out.close();
    out.open( "/tmp/F.txt" );

    for( auto&& edge : filteredEdges )
    {
      out << edge.first  << "\n"
          << edge.second << "\n\n";
    }
  }

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
