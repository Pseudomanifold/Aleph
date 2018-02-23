#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>
#include <aleph/topology/io/BipartiteAdjacencyMatrix.hh>

#include <iostream>
#include <vector>

int main( int argc, char** argv )
{
  using DataType          = double;
  using VertexType        = unsigned short;
  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  std::vector<SimplicialComplex> simplicialComplexes;
  simplicialComplexes.reserve( static_cast<unsigned>( argc - 1 ) );

  {
    aleph::topology::io::BipartiteAdjacencyMatrixReader reader;

    for( int i = 1; i < argc; i++ )
    {
      auto filename = argv[i];

      std::cerr << "* Processing " << filename << "...";

      SimplicialComplex K;
      reader( filename, K );

      std::cerr << "finished\n";

      simplicialComplexes.emplace_back( K );
    }
  }
}
