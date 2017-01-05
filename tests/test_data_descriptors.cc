#include "config/Base.hh"

#include "distances/Euclidean.hh"

#include "containers/PointCloud.hh"
#include "containers/DataDescriptors.hh"

#include "tests/Base.hh"

#include <iostream>
#include <string>
#include <vector>

using namespace aleph;

template <class D> void eccentricityTest()
{
  ALEPH_TEST_BEGIN( "Eccentricity test" );

  auto pc = load<double>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_comma_separated.txt" ) );
  auto e  = eccentricities<D>( pc );

  ALEPH_ASSERT_THROW( e.empty() == false );
  ALEPH_ASSERT_THROW( e.size()  == pc.size() );

  ALEPH_TEST_END();
}

int main()
{
  std::cerr << "-- Euclidean distance\n";

  eccentricityTest<distances::Euclidean<double> >();
}
