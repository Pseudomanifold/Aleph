/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  Given a set of persistence diagrams, it calculates their persistence
  indicator functions and---optionally---their mean indicator function
  if specified by the client.

  If not specified otherwise, all files will be written to '/tmp', and
  will have a prefix of 'PIF_'.
*/

#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <getopt.h>

#include <aleph/math/KahanSummation.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/PersistenceIndicatorFunction.hh>

#include <aleph/persistenceDiagrams/io/Raw.hh>

#include <aleph/utilities/Filesystem.hh>

using DataType           = double;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

void usage()
{
  std::cerr << "Usage: persistence_indicator_function [--mean] [--output=OUT] [--prefix=PRE] FILES\n"
            << "\n"
            << "Calculates persistence indicator functions from a set of persistence\n"
            << "diagrams, stored in FILES. Output will be written to '/tmp' and will\n"
            << "have a prefix of 'PIF_', along with the basename of the input file.\n"
            << "\n"
            << "Optionally, the mean indicator function is calculated as well, along\n"
            << "with information about the sample variance.\n"
            << "\n"
            << "Flags:\n"
            << "  -m: calculate mean persistence diagram\n"
            << "\n";
}

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "mean"  , no_argument      , nullptr, 'm' },
    { "output", required_argument, nullptr, 'o' },
    { "prefix", required_argument, nullptr, 'p' },
    { nullptr , 0                , nullptr,  0  }
  };

  if( ( argc - optind ) < 1 )
  {
    usage();
    return -1;
  }

  std::string outputDirectory = "/tmp";
  std::string prefix          = "PIF_";
  bool calculateMean = false;

  {
    int option = 0;
    while( ( option = getopt_long( argc, argv, "mo:p:", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'm':
        calculateMean = true;
        break;
      case 'o':
        outputDirectory = optarg;
        break;
      case 'p':
        prefix = optarg;
        break;
      default:
        break;
      }
    }
  }

  if( outputDirectory.empty() )
  {
    std::cerr << "* Resetting output directory to temporary directory\n";
    outputDirectory = "/tmp/";
  }

  // Ensures that the output directory ends with a slash in order to
  // indicate a proper path.
  if( outputDirectory.back() != '/' )
    outputDirectory.push_back( '/' );

  // Get filenames -----------------------------------------------------

  std::vector<std::string> filenames;
  filenames.reserve( std::size_t( argc - optind ) );

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

  using PersistenceIndicatorFunction = decltype( aleph::persistenceIndicatorFunction( PersistenceDiagram() ) );

  std::vector<PersistenceIndicatorFunction> persistenceIndicatorFunctions;
  persistenceIndicatorFunctions.reserve( persistenceDiagrams.size() );

  PersistenceIndicatorFunction mean;

  unsigned i = 0;
  for( auto&& D : persistenceDiagrams )
  {
    auto f        = aleph::persistenceIndicatorFunction( D );
    auto filename = filenames.at(i);

    if( calculateMean )
      mean += f;

    using namespace aleph::utilities;

    auto outputFilename = outputDirectory
                        + prefix
                        + stem( basename( filename ) ) + ".txt";

    std::cerr << "* Writing persistence indicator function to '" << outputFilename << "'...\n";

    std::ofstream out( outputFilename );
    out << f << "\n";

    persistenceIndicatorFunctions.emplace_back( f );

    ++i;
  }

  if( calculateMean )
  {
    mean                /= static_cast<DataType>( persistenceDiagrams.size() );
    auto outputFilename  = outputDirectory
                         + prefix
                         + "mean.txt";

    std::cerr << "* Writing mean persistence indicator function to '" << outputFilename << "'...\n";

    std::ofstream out( outputFilename );
    out << mean << "\n";

    // Since Y is supposed to be a random variable at this point, this
    // nomenclature makes sense.
    auto Y = mean.integral();

    std::cerr << "* Norm of the mean persistence indicator function: " << Y << "\n";

    std::vector<double> squaredDifferences( persistenceDiagrams.size() );

    std::transform( persistenceIndicatorFunctions.begin(), persistenceIndicatorFunctions.end(),
                    squaredDifferences.begin(),
                      [&Y] ( const PersistenceIndicatorFunction& f )
                      {
                        auto Z = f.integral();
                        return (Z-Y) * (Z-Y);
                      }
                    );

    auto s
      = aleph::math::accumulate_kahan_sorted( squaredDifferences.begin(),
                                              squaredDifferences.end(),
                                              0.0 );

    if( persistenceDiagrams.size() > 1 )
      s /= static_cast<double>( persistenceDiagrams.size() - 1 );
    else
      s = std::numeric_limits<double>::infinity();

    std::cerr << "* Sample variance: " << s << "\n";
  }
}
