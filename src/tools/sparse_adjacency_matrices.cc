#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/io/SparseAdjacencyMatrix.hh>

#include <iostream>
#include <string>
#include <vector>

int main( int argc, char** argv )
{
  using DataType          = unsigned;
  using VertexType        = std::size_t;
  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  if( argc <= 1 )
    return -1;

  std::string filename = argv[1];

  std::vector<SimplicialComplex> simplicialComplexes;
  std::vector<std::string> labels;

  aleph::topology::io::SparseAdjacencyMatrixReader reader;
  reader.setReadGraphLabels();

  std::cerr << "* Reading '" << filename << "'...";

  reader( filename, simplicialComplexes );

  std::cerr << "finished\n";
}
