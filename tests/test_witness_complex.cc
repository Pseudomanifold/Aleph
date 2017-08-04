#include <tests/Base.hh>

#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <aleph/geometry/SphereSampling.hh>
#include <aleph/geometry/WitnessComplex.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <algorithm>
#include <iterator>
#include <set>
#include <vector>

template <class SimplicialComplex> std::vector<std::size_t> bettiNumbers( SimplicialComplex K )
{
  using Simplex  = typename SimplicialComplex::ValueType;
  using DataType = typename Simplex::DataType;

  K.sort( aleph::topology::filtrations::Data<typename SimplicialComplex::ValueType>() );
  auto diagrams = aleph::calculatePersistenceDiagrams( K );

  std::vector<std::size_t> betti( diagrams.size() );

  std::transform( diagrams.begin(), diagrams.end(), betti.begin(),
                  [] ( const aleph::PersistenceDiagram<DataType>& diagram )
                  {
                    return diagram.betti();
                  } );

  return betti;
}

template <class T> void test()
{
  using Distance   = aleph::distances::Euclidean<T>;
  using PointCloud = aleph::containers::PointCloud<T>;

  PointCloud pc( 8, 2 );

  pc.set(0, {-1.0, 0.0} );
  pc.set(1, { 0.0,-1.0} );
  pc.set(2, { 1.0, 0.0} );
  pc.set(3, { 2.0, 1.0} );
  pc.set(4, { 1.0, 1.0} );
  pc.set(5, { 0.0, 2.0} );
  pc.set(6, {-1.0, 1.0} );
  pc.set(7, {-2.0, 1.0} );

  std::vector<std::size_t> indices = {0,2,4,6};

  auto K
    = aleph::geometry::buildWitnessComplex<Distance>(
        pc, indices.begin(), indices.end() );

  using Simplex    = typename decltype(K)::ValueType;
  using VertexType = typename Simplex::VertexType;

  {
    std::set<VertexType> vertices;

    K.vertices( std::inserter( vertices, vertices.begin() ) );

    ALEPH_ASSERT_EQUAL( vertices.size(), indices.size() );
  }

  auto numEdges = std::count_if( K.begin(), K.end(), [] ( const Simplex& s ) { return s.dimension() == 1; } );

  ALEPH_ASSERT_EQUAL( numEdges, 4 );
}

template <class T> void testSphereReconstruction()
{
  ALEPH_TEST_BEGIN( "Witness complexes: sphere reconstruction" );

  auto samples = aleph::geometry::sphereSampling<T>( 500 );
  auto pc      = aleph::geometry::makeSphere( samples, T(1) );

  unsigned trials = 100;
  unsigned hits   = 0;

  for( unsigned i = 0; i < trials; i++ )
  {
    std::vector<std::size_t> indices;
    aleph::geometry::generateRandomLandmarks( pc.size(), decltype(pc.size())( 12 ), std::back_inserter( indices ) );

    using Distance = aleph::distances::Euclidean<T>;

    auto K
      = aleph::geometry::buildWitnessComplex<Distance>(
          pc, indices.begin(), indices.end() );

    auto betti = bettiNumbers(K);

    if( betti == std::vector<std::size_t>( {1,0,1} ) )
      ++hits;
  }

  // TODO:
  //  - count number of hits
  //  - check that reported percentage by de Silva is reached

  for( unsigned i = 0; i < trials; i++ )
  {
    using Distance = aleph::distances::Euclidean<T>;

    std::vector<std::size_t> indices;
    aleph::geometry::generateMaxMinLandmarks( pc, 12, std::back_inserter( indices ), Distance() );

    auto K
      = aleph::geometry::buildWitnessComplex<Distance>(
          pc, indices.begin(), indices.end() );

    auto betti = bettiNumbers(K);

    if( betti == std::vector<std::size_t>( {1,0,1} ) )
      ++hits;
  }



  ALEPH_TEST_END();
}

int main(int, char**)
{
  test<float> ();
  test<double>();

  testSphereReconstruction<float> ();
  testSphereReconstruction<double>();
}
