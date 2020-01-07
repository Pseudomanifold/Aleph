#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/io/AdjacencyMatrix.hh>

#include <tests/Base.hh>

#include <sstream>

auto inputSimple =
  "0 1 3\n"
  "1 7 5\n"
  "3 5 9\n";

auto inputBroken =
  "0 1 2\n"
  "1 4\n";

template <class T> void testSimple()
{
  ALEPH_TEST_BEGIN( "Simple adjacency matrices" );

  std::stringstream stream;

  using DataType          = T;
  using Simplex           = aleph::topology::Simplex<DataType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  aleph::topology::io::AdjacencyMatrixReader reader;

  // Simple ------------------------------------------------------------

  SimplicialComplex K;

  stream.str( inputSimple );
  reader( stream, K );

  ALEPH_ASSERT_THROW( K.empty() == false );
  ALEPH_ASSERT_EQUAL( K.size(), 3 + 3 );    // 3 nodes, 3 edges

  // Ensures that the indices of the simplicial complex are consistent
  // with the dimension of the matrix.
  for( auto&& s : K )
    for( auto&& v : s )
      ALEPH_ASSERT_THROW( v <= 2 );

  std::vector<DataType> weights;

  std::transform( K.begin(), K.end(), std::back_inserter( weights ),
    [] ( const Simplex& s )
    {
      return s.data();
    }
  );

  std::sort( weights.begin(), weights.end() );

  std::vector<DataType> expectedWeights = {
    0,0,0, // vertices
    1,3,5  // edges (sorted)
  };

  ALEPH_ASSERT_THROW( weights == expectedWeights );

  // Broken ------------------------------------------------------------

  SimplicialComplex L;

  stream.clear();
  stream.str( inputBroken );

  try
  {
    reader( stream, L );
  }
  catch( std::runtime_error& )
  {
  }

  ALEPH_ASSERT_THROW( L.empty() );

  ALEPH_TEST_END();
}

int main( int, char** )
{
  testSimple<float> ();
  testSimple<double>();
}
