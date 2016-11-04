#include "config/Base.hh"

#include "containers/PointCloud.hh"

#include <iostream>
#include <string>
#include <vector>

using namespace aleph;

template <class T> void testFormats()
{
  using PointCloud = PointCloud<T>;

  std::vector<PointCloud> pointClouds;
  pointClouds.reserve( 4 );

  pointClouds.emplace_back( load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_colon_separated.txt" ) ) );
  pointClouds.emplace_back( load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_comma_separated.txt" ) ) );
  pointClouds.emplace_back( load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_space_separated.txt" ) ) );
  pointClouds.emplace_back( load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_tab_separated.txt"   ) ) );

  for( auto&& pc : pointClouds )
    std::cerr << pc.size() << "\t" << pc.dimension() << "\n";
}

int main()
{
  testFormats<double>();
  testFormats<float>();
}
