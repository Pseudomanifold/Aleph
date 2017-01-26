#include <fstream>
#include <iostream>
#include <string>

#include "filtrations/Data.hh"

#include "geometry/RipsExpander.hh"

#include "persistentHomology/ConnectedComponents.hh"

#include "topology/CliqueGraph.hh"
#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include "topology/io/EdgeLists.hh"

using DataType           = double;
using VertexType         = unsigned;
using Simplex            = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex  = aleph::topology::SimplicialComplex<Simplex>;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

int main( int argc, char** argv )
{
  if( argc <= 2 )
  {
    // TODO: Show usage
    return -1;
  }

  std::string filename = argv[1];
  unsigned k           = static_cast<unsigned>( std::stoul( argv[2] ) );

  aleph::io::EdgeListReader reader;
  reader.setReadWeights( true );
  reader.setTrimLines( true );

  SimplicialComplex K;

  {
    std::cerr << "* Reading '" << filename << "'...";

    std::ifstream in( filename );
    reader( in, K );

    std::cerr << "finished\n";
  }

  aleph::geometry::RipsExpander<SimplicialComplex> ripsExpander;
  K = ripsExpander( K, k );
  K = ripsExpander.assignMaximumWeight( K );

  K.sort( aleph::filtrations::Data<Simplex>() );

  std::cerr << "* Extracting " << k << "-cliques graph...";

  auto C
      = aleph::topology::getCliqueGraph( K, k );

  C.sort( aleph::filtrations::Data<Simplex>() );

  std::cerr << "finished\n";

  std::cerr << "* " << k << "-cliques graph has " << C.size() << " simplices\n";

  auto pd
      = aleph::calculateZeroDimensionalPersistenceDiagram( K );

  std::cout << pd << "\n";
}
