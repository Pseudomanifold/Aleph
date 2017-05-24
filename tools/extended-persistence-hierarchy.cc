#include "filtrations/Data.hh"

#include "persistenceDiagrams/PersistenceDiagram.hh"

#include "persistentHomology/ExtendedPersistenceHierarchy.hh"
#include "persistentHomology/PersistencePairing.hh"

#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include "topology/io/Function.hh"
#include "topology/io/VTK.hh"

#include "utilities/Filesystem.hh"

#include <iostream>
#include <string>
#include <vector>

// TODO: Replace this as soon as possible with a more modern option
// parser interface.
#include <getopt.h>

using DataType           = double;
using VertexType         = unsigned;
using Simplex            = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex  = aleph::topology::SimplicialComplex<Simplex>;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

void usage()
{
  // TODO: Not yet implemented
}

int main( int argc, char** argv )
{
  if( argc <= 1 )
  {
    usage();
    return -1;
  }

  std::vector<std::string> filenames;
  filenames.reserve( argc - 1 );

  for( int i = 1; i < argc; i++ )
    filenames.push_back( argv[i] );

  std::vector<SimplicialComplex> simplicialComplexes;
  simplicialComplexes.reserve( filenames.size() );

  for( auto&& filename : filenames )
  {
    std::cerr << "* Reading '" << filename << "'...";

    if( aleph::utilities::extension( filename ) == ".vtk" )
    {
      SimplicialComplex K;

      aleph::topology::io::VTKStructuredGridReader reader;
      reader( filename, K );

      K.sort( aleph::filtrations::Data<Simplex>() );

      simplicialComplexes.emplace_back( K );
    }
    else
    {
      auto complexes = aleph::topology::io::loadFunctions<SimplicialComplex>( filename );

      for( auto&& K : complexes )
        K.sort( aleph::filtrations::Data<Simplex>() );

      simplicialComplexes.insert( simplicialComplexes.end(),
                                  complexes.begin(), complexes.end() );
    }

    std::cerr << "finished\n";
  }

  for( auto&& K : simplicialComplexes )
  {
    aleph::ExtendedPersistenceHierarchy<Simplex> eph;
    auto ppe                = eph( K );   // persistence pairing & edges
    auto persistencePairing = ppe.first;  // persistence pairing
    auto edges              = ppe.second; // edges

    // TODO:
    //  - Improve output
    //  - Calculate tree-based matching
    //  - Calculate ranks

    for( auto&& edge : edges )
      std::cerr << edge.first << "--" << edge.second << "\n";
  }
}
