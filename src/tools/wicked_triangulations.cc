#include <aleph/topology/io/LexicographicTriangulation.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <iostream>
#include <string>
#include <vector>

using DataType          = bool;
using VertexType        = unsigned short;

using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

int main(int argc, char* argv[])
{
  if( argc <= 1 )
    return -1;

  std::string filename = argv[1];
  std::vector<SimplicialComplex> simplicialComplexes;

  aleph::topology::io::LexicographicTriangulationReader reader;
  reader( filename, simplicialComplexes );

  // Create missing faces ----------------------------------------------
  //
  // The triangulations are only specified by their top-level simplices,
  // so they need to be converted before being valid inputs for homology
  // calculations.

  for( auto&& K : simplicialComplexes )
  {
    K.createMissingFaces();
    K.sort();
  }

  // Calculate homology ------------------------------------------------
  //
  // We are only interested in the Betti numbers of the diagrams here as
  // the triangulations are not endowed with any weights or values.

  for( auto&& K : simplicialComplexes )
  {
    bool dualize                    = true;
    bool includeAllUnpairedCreators = true;

    auto diagrams
      = aleph::calculatePersistenceDiagrams( K,
                                            dualize,
                                            includeAllUnpairedCreators );

    for( auto&& D : diagrams )
      std::cout << D.betti() << " ";

    std::cout << "\n";
  }
}
