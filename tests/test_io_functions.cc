#include "config/Base.hh"

#include "persistentHomology/Calculation.hh"
#include "persistentHomology/ExtendedPersistenceHierarchy.hh"
#include "persistenceDiagrams/PersistenceDiagram.hh"

#include "tests/Base.hh"

#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include "topology/filtrations/Data.hh"

#include "topology/io/Function.hh"

#include <set>

template <class SimplicialComplex> aleph::PersistenceDiagram<typename SimplicialComplex::ValueType::DataType> calculatePersistenceDiagram( const SimplicialComplex& K )
{
  auto diagrams = aleph::calculatePersistenceDiagrams( K );

  ALEPH_ASSERT_EQUAL( diagrams.size(), 1 );

  return diagrams.front();
}

template <class D, class V> void test( const std::string& filename )
{
  ALEPH_TEST_BEGIN( "Functions file parsing" );

  using Simplex           = aleph::topology::Simplex<D, V>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  auto complexes
    = aleph::topology::io::loadFunctions<SimplicialComplex>( filename );

  ALEPH_ASSERT_EQUAL( complexes.size(), 2 );

  auto K = complexes.at(0);
  auto L = complexes.at(1);

  using SublevelSetFiltration   = aleph::filtrations::Data<Simplex>;
  using SuperlevelSetFiltration = aleph::filtrations::Data<Simplex, std::greater<D> >;

  ALEPH_ASSERT_EQUAL( K.size(), L.size() );
  ALEPH_ASSERT_THROW( K == L ); // modulo weights, both complexes should contain
                                // the same simplices

  K.sort( SublevelSetFiltration() );
  L.sort( SublevelSetFiltration() );

  // After sorting, the complexes must be in a different order, even
  // though their persistence pairs coincide.
  ALEPH_ASSERT_THROW( K != L );

  if( K.size() <= 9 )
  {
    auto D1 = calculatePersistenceDiagram( K );
    auto D2 = calculatePersistenceDiagram( L );

    ALEPH_ASSERT_THROW( D1 == D2 );

    aleph::ExtendedPersistenceHierarchy<Simplex> eph;
    auto edgesK = eph(K).second;
    auto edgesL = eph(L).second;

    ALEPH_ASSERT_THROW( edgesK != edgesL );
  }

  complexes
    = aleph::topology::io::loadFunctions<SimplicialComplex>(
      filename,
      [] ( D x, D y )
      {
        return std::min(x,y);
      }
  );

  ALEPH_ASSERT_EQUAL( complexes.size(), 2 );

  K = complexes.at(0);
  L = complexes.at(1);

  ALEPH_ASSERT_EQUAL( K.size(), L.size() );
  ALEPH_ASSERT_THROW( K == L ); // modulo weights, both complexes should contain
                                // the same simplices

  if( K.size() > 9 )
  {
    K.sort( SuperlevelSetFiltration() );
    L.sort( SuperlevelSetFiltration() );

    auto D1 = calculatePersistenceDiagram( K );
    auto D2 = calculatePersistenceDiagram( L );

    ALEPH_ASSERT_THROW( D1 == D2 );

    // After sorting, the complexes must be in a different order, even
    // though their persistence pairs coincide.
    ALEPH_ASSERT_THROW( K != L );

    aleph::ExtendedPersistenceHierarchy<Simplex> eph;
    auto edgesK = eph(K).second;
    auto edgesL = eph(L).second;

    ALEPH_ASSERT_THROW( edgesK != edgesL );
  }

  ALEPH_TEST_END();
}

int main()
{
  std::vector<std::string> inputs = {
    CMAKE_SOURCE_DIR + std::string( "/tests/input/Functions_simple.txt" ),
    CMAKE_SOURCE_DIR + std::string( "/tests/input/Functions_Reeb.txt" )
  };

  for( auto&& input : inputs )
  {
    test<double,unsigned>      ( input );
    test<double,unsigned short>( input );
    test<float, unsigned>      ( input );
    test<float, unsigned short>( input );
  }
}
