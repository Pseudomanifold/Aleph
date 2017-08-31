#include <aleph/config/Base.hh>

#include <tests/Base.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/io/SparseAdjacencyMatrix.hh>

#include <iostream>
#include <vector>

template <class T> void test()
{
  using DataType          = float;
  using VertexType        = T;
  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  std::vector<SimplicialComplex> complexes;

  aleph::topology::io::SparseAdjacencyMatrixReader reader;
  reader.setReadEdgeAttributes();
  reader.setEdgeAttributeIndex(0);
  reader( CMAKE_SOURCE_DIR + std::string( "/tests/input/Simple_adjacency_matrix_A.txt"), complexes );

  ALEPH_ASSERT_THROW( complexes.empty() == false );
  ALEPH_ASSERT_EQUAL( complexes.size(), 3 );
  ALEPH_ASSERT_EQUAL( complexes[0].size(), 6 );
  ALEPH_ASSERT_EQUAL( complexes[1].size(), 3 );
  ALEPH_ASSERT_EQUAL( complexes[2].size(), 3 );
}

int main(int, char**)
{
  test<unsigned>();
}
