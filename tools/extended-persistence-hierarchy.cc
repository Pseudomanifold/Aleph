
#include "persistenceDiagrams/PersistenceDiagram.hh"

#include "persistentHomology/ExtendedPersistenceHierarchy.hh"
#include "persistentHomology/PersistencePairing.hh"

#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include "topology/filtrations/Data.hh"

#include "topology/io/Function.hh"
#include "topology/io/VTK.hh"

#include "utilities/Filesystem.hh"

#include <functional>
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
  std::cerr << "Usage: extended-persistence-hierarchy [--superlevels] [--sublevels] FILES\n";
}

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "superlevels", no_argument, nullptr, 'S' },
    { "sublevels"  , no_argument, nullptr, 's' },
    { nullptr      , 0          , nullptr,  0  }
  };

  bool calculateSuperlevelSets = false;

  int option = 0;
  while( ( option = getopt_long( argc, argv, "Ss", commandLineOptions, nullptr ) ) != -1 )
  {
    switch( option )
    {
    case 'S':
      calculateSuperlevelSets = true;
      break;
    case 's':
      calculateSuperlevelSets = false;
      break;
    default:
      break;
    }
  }

  if( ( argc - optind ) <= 0 )
  {
    usage();
    return -1;
  }

  std::vector<std::string> filenames;
  filenames.reserve( argc - optind );

  for( int i = optind; i < argc; i++ )
    filenames.push_back( argv[i] );

  std::vector<SimplicialComplex> simplicialComplexes;
  simplicialComplexes.reserve( filenames.size() );

  // Select a functor for calculating the weights when reading
  // a simplicial complex.
  auto functor = calculateSuperlevelSets
               ? [] ( DataType a, DataType b ) { return std::min(a,b); }
               : [] ( DataType a, DataType b ) { return std::max(a,b); };

  for( auto&& filename : filenames )
  {
    std::cerr << "* Reading '" << filename << "'...";

    if( aleph::utilities::extension( filename ) == ".vtk" )
    {
      SimplicialComplex K;

      aleph::topology::io::VTKStructuredGridReader reader;
      reader( filename, K, functor );

      if( calculateSuperlevelSets )
        K.sort( aleph::filtrations::Data<Simplex, std::greater<DataType> >() );
      else
        K.sort( aleph::filtrations::Data<Simplex, std::less<DataType> >() );

      simplicialComplexes.emplace_back( K );
    }
    else
    {
      auto complexes
        = aleph::topology::io::loadFunctions<SimplicialComplex>( filename,
                                                                 functor );

      for( auto&& K : complexes )
      {
        if( calculateSuperlevelSets )
          K.sort( aleph::filtrations::Data<Simplex, std::less<DataType> >() );
        else
          K.sort( aleph::filtrations::Data<Simplex, std::greater<DataType> >() );
      }

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
