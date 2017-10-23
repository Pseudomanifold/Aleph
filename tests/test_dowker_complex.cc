#include <tests/Base.hh>

#include <aleph/geometry/DowkerComplex.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <vector>

template <class T> void test()
{
  ALEPH_TEST_BEGIN( "Simple directed networks" );

  using Matrix = std::vector< std::vector<T> >;

  Matrix X( 3, std::vector<T>( 3 ) );

  X[0][1] = T(6);
  X[0][2] = T(4);
  X[1][0] = T(1);
  X[1][2] = T(5);
  X[2][0] = T(3);
  X[2][1] = T(3);

  Matrix Y( 3, std::vector<T>( 3 ) );

  Y[0][1] = T(6);
  Y[0][2] = T(3);
  Y[1][0] = T(1);
  Y[1][2] = T(5);
  Y[2][0] = T(4);
  Y[2][1] = T(3);

  auto X_pairs = aleph::geometry::admissiblePairs( X, T(6) );
  auto Y_pairs = aleph::geometry::admissiblePairs( Y, T(6) );

  ALEPH_ASSERT_THROW( X_pairs.empty() == false );
  ALEPH_ASSERT_THROW( Y_pairs.empty() == false );

  ALEPH_ASSERT_EQUAL( X.size(), Y.size() );

  auto X_dowkerComplexes = aleph::geometry::buildDowkerSinkSourceComplexes<unsigned, T>( X_pairs );
  auto Y_dowkerComplexes = aleph::geometry::buildDowkerSinkSourceComplexes<unsigned, T>( Y_pairs );

  auto X_source = X_dowkerComplexes.first;
  auto X_sink   = X_dowkerComplexes.second;
  auto Y_source = Y_dowkerComplexes.first;
  auto Y_sink   = Y_dowkerComplexes.second;

  ALEPH_ASSERT_EQUAL( X_source.size(), Y_source.size() );
  ALEPH_ASSERT_EQUAL( X_sink.size(),   Y_sink.size()   );

  auto X_source_diagrams = aleph::calculatePersistenceDiagrams( X_source );
  auto X_sink_diagrams   = aleph::calculatePersistenceDiagrams( X_sink   );
  auto Y_source_diagrams = aleph::calculatePersistenceDiagrams( Y_source );
  auto Y_sink_diagrams   = aleph::calculatePersistenceDiagrams( Y_sink   );

  ALEPH_ASSERT_EQUAL( X_source_diagrams.size(), X_sink_diagrams.size() );
  ALEPH_ASSERT_EQUAL( Y_source_diagrams.size(), Y_sink_diagrams.size() );

  for( auto&& D : X_source_diagrams )
    D.removeDiagonal();

  for( auto&& D : X_sink_diagrams )
    D.removeDiagonal();

  for( auto&& D : Y_source_diagrams )
    D.removeDiagonal();

  for( auto&& D : Y_sink_diagrams )
    D.removeDiagonal();

  ALEPH_ASSERT_EQUAL( X_source_diagrams.size(), Y_source_diagrams.size() );
  ALEPH_ASSERT_EQUAL( X_sink_diagrams.size()  , Y_sink_diagrams.size() );

  ALEPH_ASSERT_EQUAL( X_source_diagrams.size(), 2 );
  ALEPH_ASSERT_EQUAL( Y_source_diagrams.size(), 2 );

  ALEPH_ASSERT_THROW( X_source_diagrams.front() == Y_source_diagrams.front() );
  ALEPH_ASSERT_THROW( X_source_diagrams.back()  != Y_source_diagrams.back() );

  ALEPH_TEST_END();
}

int main( int, char** )
{
  test<float> ();
  test<double>();
}
