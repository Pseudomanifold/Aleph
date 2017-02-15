#include <algorithm>
#include <fstream>
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

int main( int argc, char** argv )
{
  if( argc <= 2 )
  {
    // TODO: Show usage
    return -1;
  }

  std::string filename = argv[1];
  double threshold     = std::stod( argv[2] );

  SimplicialComplex K;

  // Input -------------------------------------------------------------
  //
  // TODO: This is copied from the clique persistence diagram
  // calculation. It would make sense to share some functions
  // between the two applications.

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

  // Thresholding ------------------------------------------------------

  {
    std::cerr << "* Filtering input data to threshold epsilon=" << threshold << "...";

    std::vector<Simplex> simplices;

    std::remove_copy_if( K.begin(), K.end(), std::back_inserter( simplices ),
                         [&threshold] ( const Simplex& s )
                         {
                           return s.data() > threshold;
                         } );

    K = SimplicialComplex( simplices.begin(), simplices.end() );

    std::cerr << "finished\n";
  }

  // Expansion ---------------------------------------------------------

  unsigned maxK = 12;

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

    // TODO:
    //  - Calculate connected components of clique graph
    //  - Look up original simplices in the simplicial complex
    //  - Create output
  }
}
