#include <aleph/persistenceDiagrams/Norms.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistentHomology/Calculation.hh>

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

  // 1. Read simplicial complexes --------------------------------------

  std::vector<SimplicialComplex> simplicialComplexes;
  simplicialComplexes.reserve( static_cast<unsigned>( argc - 1 ) );

  {
    aleph::topology::io::BipartiteAdjacencyMatrixReader reader;

    for( int i = 1; i < argc; i++ )
    {
      auto filename = std::string( argv[i] );

      std::cerr << "* Processing " << filename << "...";

      SimplicialComplex K;
      reader( filename, K );

      std::cerr << "finished\n";

      simplicialComplexes.emplace_back( K );
    }
  }

  // 2. Calculate persistent homology ----------------------------------

  for( std::size_t i = 0; i < simplicialComplexes.size(); i++ )
  {
    auto&& K      = simplicialComplexes[i];
    auto diagrams = aleph::calculatePersistenceDiagrams( K );
    auto&& D      = diagrams.front();

    D.removeUnpaired();

    std::cout << i << "\t" << aleph::pNorm( D ) << "\n";
  }
}
