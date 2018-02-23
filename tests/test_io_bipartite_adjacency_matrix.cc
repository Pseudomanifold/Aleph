#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/io/BipartiteAdjacencyMatrix.hh>

#include <tests/Base.hh>

#include <sstream>

auto inputSimple =
  "0 1 2\n"
  "3 4 5\n";

auto inputBroken =
  "0 1 2\n"
  "1 4\n";

auto inputHorror =
  "1 1 1\n"
  "1 1\n";

template <class T> void testSimple()
{
  ALEPH_TEST_BEGIN( "Simple adjacency matrices" );

  std::stringstream stream;

  using DataType          = T;
  using Simplex           = aleph::topology::Simplex<DataType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  aleph::topology::io::BipartiteAdjacencyMatrixReader reader;

  // Simple ------------------------------------------------------------

  SimplicialComplex K;

  stream.str( inputSimple );
  reader( stream, K );

  ALEPH_ASSERT_THROW( K.empty() == false );
  ALEPH_ASSERT_EQUAL( K.size(), 6 + 5 );    // 6 nodes, 5 edges

  std::vector<DataType> weights;

  std::transform( K.begin(), K.end(), std::back_inserter( weights ),
    [] ( const Simplex& s )
    {
      return s.data();
    }
  );

  std::sort( weights.begin(), weights.end() );

  std::vector<DataType> expectedWeights = {
    0,0,0,0,0,0,  // nodes
    1,2,3,4,5     // edges (sorted)
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

  // Horror ------------------------------------------------------------

  SimplicialComplex M;

  stream.clear();
  stream.str( inputHorror );

  try
  {
    reader( stream, M );
  }
  catch( std::runtime_error& )
  {
  }

  ALEPH_ASSERT_THROW( M.empty() );

  ALEPH_TEST_END();
}

int main( int, char** )
{
  testSimple<float> ();
  testSimple<double>();
}
