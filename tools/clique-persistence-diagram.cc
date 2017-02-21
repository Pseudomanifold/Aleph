#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>

#include <cmath>

#include "filtrations/Data.hh"

#include "geometry/RipsExpander.hh"

#include "persistentHomology/ConnectedComponents.hh"

#include "topology/CliqueGraph.hh"
#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include "topology/io/EdgeLists.hh"
#include "topology/io/GML.hh"

#include "utilities/Filesystem.hh"

using DataType           = double;
using VertexType         = unsigned;
using Simplex            = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex  = aleph::topology::SimplicialComplex<Simplex>;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

std::string formatOutput( const std::string& prefix, unsigned k, unsigned K )
{
  std::ostringstream stream;
  stream << prefix;
  stream << std::setw( int( std::log10( K ) + 1 ) ) << std::setfill( '0' ) << k;
  stream << ".txt";

  return stream.str();
}

int main( int argc, char** argv )
{
  if( argc <= 2 )
  {
    // TODO: Show usage
    return -1;
  }

  std::string filename = argv[1];
  unsigned maxK        = static_cast<unsigned>( std::stoul( argv[2] ) );

  SimplicialComplex K;

  // Input -------------------------------------------------------------

  std::cerr << "* Reading '" << filename << "'...";

  if( aleph::utilities::extension( filename ) == ".gml" )
  {
    aleph::topology::io::GMLReader reader;
    reader( filename, K );
  }
  else
  {
    aleph::io::EdgeListReader reader;
    reader.setReadWeights( true );
    reader.setTrimLines( true );

    reader( filename, K );
  }

  std::cerr << "finished\n";

  // Expansion ---------------------------------------------------------

  aleph::geometry::RipsExpander<SimplicialComplex> ripsExpander;
  K = ripsExpander( K, maxK );
  K = ripsExpander.assignMaximumWeight( K );

  K.sort( aleph::filtrations::Data<Simplex>() );

  for( unsigned k = 1; k <= maxK; k++ )
  {
    std::cerr << "* Extracting " << k << "-cliques graph...";

    auto C
        = aleph::topology::getCliqueGraph( K, k );

    C.sort( aleph::filtrations::Data<Simplex>() );

    std::cerr << "finished\n";

    std::cerr << "* " << k << "-cliques graph has " << C.size() << " simplices\n";

    if( C.empty() )
    {
      std::cerr << "* Stopping here because no further cliques for processing exist\n";
      break;
    }

    auto pd
        = aleph::calculateZeroDimensionalPersistenceDiagram( C );

    {
      using namespace aleph::utilities;
      auto outputFilename = formatOutput( "/tmp/" + stem( basename( filename ) ) + "_k", k, maxK );

      std::cerr << "* Storing output in '" << outputFilename << "'...\n";

      pd.removeDiagonal();

      std::ofstream out( outputFilename );
      out << "# Original filename: " << filename << "\n";
      out << "# k                : " << k        << "\n";
      out << pd << "\n";
    }
  }
}
