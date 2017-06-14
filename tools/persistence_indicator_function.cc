/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  Given a set of persistence diagrams, it calculates their persistence
  indicator functions and---optionally---their mean indicator function
  if specified by the client.

  All files will be written to '/tmp', prefixed with 'PIF_'.
*/

#include <iostream>
#include <string>
#include <vector>

#include <getopt.h>

#include "persistenceDiagrams/PersistenceDiagram.hh"
#include "persistenceDiagrams/PersistenceIndicatorFunction.hh"

#include "persistenceDiagrams/io/Raw.hh"

#include "utilities/Filesystem.hh"

using DataType           = double;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "mean" , no_argument, nullptr, 'm' },
    { nullptr, 0          , nullptr,  0  }
  };

  if( ( argc - optind ) < 1 )
  {
    // TODO: Show usage
    return -1;
  }

  bool calculateMean = false;

  {
    int option = 0;
    while( ( option = getopt_long( argc, argv, "m", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'm':
        calculateMean = true;
        break;
      default:
        break;
      }
    }
  }

  // Get filenames -----------------------------------------------------

  std::vector<std::string> filenames;
  filenames.reserve( argc - optind );

  for( int i = optind; i < argc; i++ )
    filenames.push_back( argv[i] );

  // Load persistence diagrams -----------------------------------------

  std::vector<PersistenceDiagram> persistenceDiagrams;
  persistenceDiagrams.reserve( filenames.size() );

  for( auto&& filename : filenames )
  {
    std::cerr << "* Processing '" << filename << "'...";

    PersistenceDiagram persistenceDiagram = aleph::io::load<DataType>( filename );

    // FIXME: This is only required in order to ensure that the
    // persistence indicator function has a finite integral; it
    // can be solved more elegantly by using a special value to
    // indicate infinite intervals.
    persistenceDiagram.removeUnpaired();

    persistenceDiagrams.push_back( persistenceDiagram );

    std::cerr << "finished\n";
  }

  // Calculate persistence indicator functions -------------------------

  decltype( aleph::persistenceIndicatorFunction( PersistenceDiagram() ) ) mean;

  unsigned i = 0;
  for( auto&& D : persistenceDiagrams )
  {
    auto f        = aleph::persistenceIndicatorFunction( D );
    auto filename = filenames.at(i);

    if( calculateMean )
      mean += f;

    using namespace aleph::utilities;

    auto outputFilename = "/tmp/PIF_" + stem( basename( filename ) ) + ".txt";

    std::cerr << "* Writing persistence indicator function to '" << outputFilename << "'...\n";

    std::ofstream out( outputFilename );
    out << f << "\n";

    ++i;
  }

  if( calculateMean )
  {
    mean                /= static_cast<DataType>( persistenceDiagrams.size() );
    auto outputFilename  = "/tmp/PIF_mean.txt";

    std::cerr << "* Writing mean persistence indicator function to '" << outputFilename << "'...\n";

    std::ofstream out( outputFilename );
    out << mean << "\n";
  }
}
