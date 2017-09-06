#include <tests/Base.hh>

#include <aleph/config/Eigen.hh>

#include <aleph/geometry/HeatKernel.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

template <class T> aleph::topology::SimplicialComplex< aleph::topology::Simplex<T, unsigned> > createTestSimplicialComplex()
{
  using Simplex           = typename aleph::topology::Simplex<T, unsigned>;
  using SimplicialComplex = typename aleph::topology::SimplicialComplex<Simplex>;

  std::vector<Simplex> simplices = {
    Simplex( {1},  T(1) ),
    Simplex( {2},  T(1) ),
    Simplex( {3},  T(1) ),
    Simplex( {4},  T(1) ),
    Simplex( {1,2},T(1) ),
    Simplex( {1,3},T(1) ),
    Simplex( {2,3},T(1) ),
    Simplex( {2,4},T(1) )
  };

  return SimplicialComplex( simplices.begin(), simplices.end() );
}

#ifdef ALEPH_WITH_EIGEN

template <class T> void testWeightedLaplacianMatrix()
{
  ALEPH_TEST_BEGIN( "Weighted Laplacian matrix" );

  auto K = createTestSimplicialComplex<T>();
  auto M = aleph::geometry::weightedAdjacencyMatrix( K );

  for( auto&& i : {0,1,2,3} )
    ALEPH_ASSERT_EQUAL( M(i,i), T(0) );

  ALEPH_ASSERT_EQUAL( M(0,1), T(1) );
  ALEPH_ASSERT_EQUAL( M(0,2), T(1) );
  ALEPH_ASSERT_EQUAL( M(1,2), T(1) );
  ALEPH_ASSERT_EQUAL( M(1,3), T(1) );

  for( auto&& i : {0,1,2,3} )
    for( auto&& j : {0,1,2,3} )
      ALEPH_ASSERT_EQUAL( M(i,j), M(j,i) );

  auto L = aleph::geometry::weightedLaplacianMatrix( K );

  std::vector<T> actualValues;
  actualValues.reserve( std::size_t( L.size() ) );

  for( unsigned i = 0; i < L.rows(); i++ )
    for( unsigned j = 0; j < L.cols(); j++ )
      actualValues.push_back( L(i,j) );

  std::vector<T> expectedValues = {
      2,-1,-1, 0,
     -1, 3,-1,-1,
     -1,-1, 2, 0,
      0,-1, 0, 1,
  };

  ALEPH_ASSERT_THROW( actualValues == expectedValues );

  ALEPH_TEST_END();
}

template <class T> void testHeatKernelSimple()
{
  ALEPH_TEST_BEGIN( "Simple heat kernel test" );

  auto K = createTestSimplicialComplex<T>();
  aleph::geometry::HeatKernel hk( K );

  auto samples1 = hk.logarithmicSamplingInterval( 20 );

  for( auto&& t : samples1 )
  {
    auto trace = hk.trace(t);
    std::cerr << t << "\t" << trace << "\n";
  }

  ALEPH_TEST_END();
}

#endif

int main( int, char** )
{
#ifdef ALEPH_WITH_EIGEN
  testWeightedLaplacianMatrix<float> ();
  testWeightedLaplacianMatrix<double>();

  testHeatKernelSimple<float> ();
  testHeatKernelSimple<double>();
#endif
}
