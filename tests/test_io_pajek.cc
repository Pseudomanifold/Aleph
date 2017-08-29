#include <aleph/config/Base.hh>

#include <tests/Base.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/io/Pajek.hh>
#include <aleph/topology/io/SimplicialComplexReader.hh>

#include <set>

template <class D, class V> void test( const std::string& filename )
{
  ALEPH_TEST_BEGIN( "Pajek file parsing" );

  using Simplex           = aleph::topology::Simplex<D, V>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K;

  aleph::topology::io::PajekReader reader;
  reader( filename, K );

  ALEPH_ASSERT_EQUAL( K.size(), 22 );
  ALEPH_ASSERT_EQUAL( std::count_if( K.begin(), K.end(), [] (const Simplex& s) { return s.dimension() == 0; } ), 10 );
  ALEPH_ASSERT_EQUAL( std::count_if( K.begin(), K.end(), [] (const Simplex& s) { return s.dimension() == 1; } ), 12 );

  std::set<V> vertices;
  K.vertices( std::inserter( vertices, vertices.begin() ) );

  ALEPH_ASSERT_EQUAL( *vertices.begin(),  1  );
  ALEPH_ASSERT_EQUAL( *vertices.rbegin(), 10 );

  auto itSigma = K.find( Simplex( {2,8} ) );
  auto itTau   = K.find( Simplex( {5,7} ) );

  ALEPH_ASSERT_THROW( itSigma != K.end() );
  ALEPH_ASSERT_THROW( itTau   != K.end() );

  auto w1 = itSigma->data();
  auto w2 = itTau->data();

  if( w1 == w2 )
  {
    ALEPH_ASSERT_EQUAL( w1, 0 );
    ALEPH_ASSERT_EQUAL( w2, 0 );
  }
  else
  {
    ALEPH_ASSERT_EQUAL( w1, 23 );
    ALEPH_ASSERT_EQUAL( w2, 42 );
  }

  {
    SimplicialComplex L;

    aleph::topology::io::SimplicialComplexReader reader;
    reader( filename, L );

    ALEPH_ASSERT_THROW( K == L );

    auto sigma1 = *K.find( Simplex( {2,8} ) );
    auto sigma2 = *L.find( Simplex( {2,8} ) );

    ALEPH_ASSERT_EQUAL( sigma1.data(), sigma2.data() );
  }

  ALEPH_TEST_END();
}

int main()
{
  std::vector<std::string> inputs = {
    CMAKE_SOURCE_DIR + std::string( "/tests/input/Simple.net" ),
    CMAKE_SOURCE_DIR + std::string( "/tests/input/Simple_with_labels.net" )
  };

  for( auto&& input : inputs )
  {
    test<double,unsigned>      ( input );
    test<double,unsigned short>( input );
    test<float, unsigned>      ( input );
    test<float, unsigned short>( input );
  }
}
