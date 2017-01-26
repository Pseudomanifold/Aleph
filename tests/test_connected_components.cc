#include "config/Base.hh"

#include "containers/PointCloud.hh"

#include "distances/Euclidean.hh"

#include "filtrations/Data.hh"

#include "geometry/BruteForce.hh"
#include "geometry/RipsSkeleton.hh"

#include "tests/Base.hh"

#include "persistentHomology/Calculation.hh"

#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include <vector>

using namespace aleph::geometry;
using namespace aleph::topology;
using namespace aleph;

template <class T> void test()
{
  ALEPH_TEST_BEGIN( "Point cloud loading" );

  using PointCloud = PointCloud<T>;
  using Distance   = aleph::distances::Euclidean<T>;

  PointCloud pointCloud = load<T>( CMAKE_SOURCE_DIR + std::string( "/tests/input/Iris_colon_separated.txt" ) );

  ALEPH_ASSERT_THROW( pointCloud.size()      == 150 );
  ALEPH_ASSERT_THROW( pointCloud.dimension() ==   4);

  ALEPH_TEST_END();

  ALEPH_TEST_BEGIN( "Rips skeleton calculation" );

  using Wrapper      = BruteForce<PointCloud, Distance>;
  using RipsSkeleton = RipsSkeleton<Wrapper>;

  Wrapper wrapper( pointCloud );
  RipsSkeleton ripsSkeleton;

  auto K = ripsSkeleton( wrapper, 1.0 );


  using Simplex = typename decltype(K)::ValueType;

  ALEPH_ASSERT_THROW( K.empty() == false );
  ALEPH_ASSERT_THROW( std::count_if( K.begin(), K.end(), [] ( const Simplex& s ) { return s.dimension() == 1; } ) != 0 );

  ALEPH_TEST_END();

  ALEPH_TEST_BEGIN( "Zero-dimensional persistent homology calculation" );

  K.sort( filtrations::Data<Simplex>() );

  auto diagrams = calculatePersistenceDiagrams( K );

  ALEPH_ASSERT_THROW( diagrams.empty() == false );

  auto diagram1 = diagrams.front();

  ALEPH_ASSERT_THROW( diagram1.empty() == false );
  ALEPH_ASSERT_THROW( diagram1.size()  == pointCloud.size() );

  ALEPH_TEST_END();
}

int main()
{
  test<float> ();
  test<double>();
}
