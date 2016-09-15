#include "io/PLY.hh"
#include "utilities/Timer.hh"

#include "PersistentHomologyCalculation.hh"
#include "PersistenceDiagramNorms.hh"

#include <iostream>

using DataType   = double;
using VertexType = unsigned;

int main( int argc, char** argv )
{
  std::string filename;
  std::string property = "quality";

  if( argc == 1 )
    return -1;

  if( argc >= 2 )
    filename = argv[1];

  if( argc >= 3 )
    property = argv[2];

  auto K
    = aleph::io::loadPLY<DataType, VertexType>( filename, property );

  std::cerr << "* Loaded simplicial complex with " << K.size() << " simplices\n";

  aleph::utilities::Timer timer;

  auto diagrams
    = aleph::calculatePersistenceDiagrams( K );

  std::cerr << "* Calculated " << diagrams.size() << " persistence diagrams in " << timer.elapsed_s() << "s\n";

  for( auto&& D : diagrams )
    std::cout << D << "\n";

  for( auto&& D : diagrams )
  {
    std::cerr << "* Total degree-1 persistence: " << aleph::totalPersistence( D, 1.0 ) << "\n"
              << "* Total degree-2 persistence: " << aleph::totalPersistence( D, 2.0 ) << "\n"
              << "* 1-norm:                     " << aleph::pNorm( D, 1.0 ) << "\n"
              << "* 2-norm:                     " << aleph::pNorm( D, 2.0 ) << "\n";
  }
}
