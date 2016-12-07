#include "config/Base.hh"

#include "containers/PointCloud.hh"

#include "distances/Euclidean.hh"

#include "geometry/BruteForce.hh"
#include "geometry/RipsExpander.hh"
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

  ALEPH_TEST_BEGIN( "Vietoris--Rips expansion" );

  using Wrapper      = BruteForce<PointCloud, Distance>;
  using RipsSkeleton = RipsSkeleton<Wrapper>;

  Wrapper wrapper( pointCloud );
  RipsSkeleton ripsSkeleton;

  auto K = ripsSkeleton( wrapper, 1.0 );

  using Simplex = typename decltype(K)::ValueType;

  ALEPH_ASSERT_THROW( K.empty() == false );
  ALEPH_ASSERT_THROW( std::count_if( K.begin(), K.end(), [] ( const Simplex& s ) { return s.dimension() == 1; } ) != 0 );

  bool dualize     = true;
  bool notDualized = false;

  auto diagrams1 = calculatePersistenceDiagrams( K, dualize );
  auto diagrams2 = calculatePersistenceDiagrams( K, notDualized );

  ALEPH_ASSERT_THROW( diagrams1.size() == diagrams2.size() );

  for( std::size_t i = 0; i < diagrams1.size(); i++ )
  {
    auto&& D1 = diagrams1.at(i);
    auto&& D2 = diagrams2.at(i);

    ALEPH_ASSERT_THROW( D1.dimension() == D2.dimension() );
    ALEPH_ASSERT_THROW( D1 == D2 );
  }

  ALEPH_TEST_END();
}

int main()
{
  test<float> ();
  test<double>();
}
