#include "algorithms/Standard.hh"
#include "algorithms/Twist.hh"

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

#include "topology/representations/List.hh"
#include "topology/representations/Set.hh"
#include "topology/representations/Vector.hh"

#include <vector>

using namespace aleph::persistentHomology::algorithms;
using namespace aleph::geometry;
using namespace aleph::topology;
using namespace aleph;

template <class R, class S> auto testInternal( const S& K ) -> std::vector< PersistenceDiagram<typename S::ValueType::DataType> >
{
  using Simplex  = typename S::ValueType;
  using DataType = typename Simplex::DataType;

  std::vector< PersistenceDiagram<DataType> > diagrams;

  bool dualize     = true;
  bool notDualized = false;

  auto diagrams1 = calculatePersistenceDiagrams<Standard, R>( K, dualize );
  auto diagrams2 = calculatePersistenceDiagrams<Standard, R>( K, notDualized );
  auto diagrams3 = calculatePersistenceDiagrams<Twist, R>( K, dualize );
  auto diagrams4 = calculatePersistenceDiagrams<Twist, R>( K, notDualized );

  ALEPH_ASSERT_THROW( diagrams1.size() == diagrams2.size() );
  ALEPH_ASSERT_THROW( diagrams2.size() == diagrams3.size() );
  ALEPH_ASSERT_THROW( diagrams3.size() == diagrams4.size() );

  for( std::size_t i = 0; i < diagrams1.size(); i++ )
  {
    auto&& D1 = diagrams1.at(i);
    auto&& D2 = diagrams2.at(i);
    auto&& D3 = diagrams3.at(i);
    auto&& D4 = diagrams4.at(i);

    ALEPH_ASSERT_THROW( D1.dimension() == D2.dimension() );
    ALEPH_ASSERT_THROW( D2.dimension() == D3.dimension() );
    ALEPH_ASSERT_THROW( D3.dimension() == D4.dimension() );
    ALEPH_ASSERT_THROW( D1 == D2 );
    ALEPH_ASSERT_THROW( D2 == D3 );
    ALEPH_ASSERT_THROW( D3 == D4 );
  }

  diagrams.insert( diagrams.end(), diagrams1.begin(), diagrams1.end() );
  diagrams.insert( diagrams.end(), diagrams2.begin(), diagrams2.end() );
  diagrams.insert( diagrams.end(), diagrams3.begin(), diagrams3.end() );
  diagrams.insert( diagrams.end(), diagrams4.begin(), diagrams4.end() );

  return diagrams;
}

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
  using Index   = typename Simplex::VertexType;

  ALEPH_ASSERT_THROW( K.empty() == false );
  ALEPH_ASSERT_THROW( std::count_if( K.begin(), K.end(), [] ( const Simplex& s ) { return s.dimension() == 1; } ) != 0 );

  auto diagrams3 = testInternal<representations::List<Index> >( K );
  auto diagrams1 = testInternal<representations::Set<Index> >( K );
  auto diagrams2 = testInternal<representations::Vector<Index> >( K );

  ALEPH_ASSERT_THROW( diagrams1.size() == diagrams2.size() );
  ALEPH_ASSERT_THROW( diagrams2.size() == diagrams3.size() );

  for( std::size_t i = 0; i < diagrams1.size(); i++ )
  {
    auto&& D1 = diagrams1.at(i);
    auto&& D2 = diagrams2.at(i);
    auto&& D3 = diagrams3.at(i);

    ALEPH_ASSERT_THROW( D1.dimension() == D2.dimension() );
    ALEPH_ASSERT_THROW( D2.dimension() == D3.dimension() );
    ALEPH_ASSERT_THROW( D1 == D2 );
    ALEPH_ASSERT_THROW( D2 == D3 );
  }

  ALEPH_TEST_END();
}

int main()
{
  test<float> ();
  test<double>();
}
