#include <tests/Base.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>
#include <aleph/topology/io/LexicographicTriangulation.hh>

#include <iterator>
#include <sstream>
#include <vector>

template <class T> void test()
{
  std::stringstream stream;
  stream << "manifold_2_6_3=[  [1,2,3],[1,2,4],[1,3,5],[1,4,5],\n"
         << "[2,3,6],\n"
         << "[2,4,6],[3,5,6],[4,5,6]\n"
         << "]\n\n"
         << "manifold_simple = [ [1,2,3,] ]\n";

  using DataType          = bool;
  using VertexType        = T;
  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  std::vector<SimplicialComplex> simplicialComplexes;

  aleph::topology::io::LexicographicTriangulationReader reader;

  // ugly!
  reader.operator()<SimplicialComplex>( stream,
                                        std::back_inserter( simplicialComplexes ) );

  ALEPH_ASSERT_EQUAL( simplicialComplexes.size(), 2 );

  auto&& K = simplicialComplexes.front();
  auto&& L = simplicialComplexes.back();

  ALEPH_ASSERT_THROW( K.empty() == false );
  ALEPH_ASSERT_THROW( L.empty() == false );

  ALEPH_ASSERT_EQUAL( K.size(), 8 );
  ALEPH_ASSERT_EQUAL( L.size(), 1 );
}

int main(int, char**)
{
  test<unsigned>();
}
