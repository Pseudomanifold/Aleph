#include "config/Base.hh"

#include "containers/PointCloud.hh"

#include "tests/Base.hh"

#include <iostream>
#include <string>
#include <vector>

using namespace aleph;

template <class T> void testFormats()
{
  ALEPH_TEST_BEGIN( "Point cloud formats" );

  using PointCloud = PointCloud<T>;

  std::vector<PointCloud> pointClouds;
  pointClouds.reserve( 4 );

  pointClouds.emplace_back( load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_colon_separated.txt" ) ) );
  pointClouds.emplace_back( load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_comma_separated.txt" ) ) );
  pointClouds.emplace_back( load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_space_separated.txt" ) ) );
  pointClouds.emplace_back( load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_tab_separated.txt"   ) ) );

  for( auto&& pc : pointClouds )
  {
    ALEPH_ASSERT_THROW( pc.size()      == 150 );
    ALEPH_ASSERT_THROW( pc.dimension() == 4 );
    ALEPH_ASSERT_THROW( pc.empty()     == false );
  }

  for( auto&& pc1 : pointClouds )
    for( auto&& pc2 : pointClouds )
      ALEPH_ASSERT_THROW( pc1 == pc2 );

  ALEPH_TEST_END();
}

template <class T> void testAccess()
{
  ALEPH_TEST_BEGIN( "Point cloud access" );

  auto pc
    = load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_comma_separated.txt" ) );

  {
    std::vector<T> actual;
    std::vector<T> expected = { static_cast<T>( 5.9 ),
                                static_cast<T>( 3.0 ),
                                static_cast<T>( 5.1 ),
                                static_cast<T>( 1.8 ) };

    pc.get( 149, std::back_inserter( actual ) );

    ALEPH_ASSERT_THROW( actual == expected );
  }

  {
    std::vector<T> p = { T(1), T(2), T(3), T(4) };
    std::vector<T> q;

    pc.set( 149, p.begin(), p.end() );
    pc.get( 149, std::back_inserter( q ) );

    ALEPH_ASSERT_THROW( p == q );
  }

  {
    std::vector<T> p1 = { T(1), T(2), T(3), T(4) };
    std::vector<T> p2;
    std::vector<T> q1 = { T(5), T(6), T(7), T(8) };
    std::vector<T> q2;

    pc.set( 148, {1,2,3,4} );
    pc.set( 149, {5,6,7,8} );

    pc.get( 148, std::back_inserter( p2 ) );
    pc.get( 149, std::back_inserter( q2 ) );

    ALEPH_ASSERT_THROW( p1 == p2 );
    ALEPH_ASSERT_THROW( q1 == q2 );
  }

  {
    std::vector<T> p;

    ALEPH_EXPECT_EXCEPTION(
      pc.get( 151, std::back_inserter( p ) ),
      std::runtime_error
    );
  }

  {
    std::vector<T> p;

    ALEPH_EXPECT_EXCEPTION(
      pc.set( 151, {1,2,3} ),
      std::runtime_error
    );
  }

  ALEPH_TEST_END();
}

int main()
{
  std::cerr << "-- float\n";

  testFormats<float> ();
  testAccess<float>  ();

  std::cerr << "-- double\n";

  testFormats<double>();
  testAccess<double> ();
}
