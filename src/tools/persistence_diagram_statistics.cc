/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  It analyses various aspects, such as the total persistence, of a set
  of persistence diagrams and writes all statistics to STDOUT. Results
  are formatted as comma-separated values (CSV).

  Original author: Bastian Rieck
*/

#include <aleph/persistenceDiagrams/Extraction.hh>
#include <aleph/persistenceDiagrams/Norms.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistenceDiagrams/io/Raw.hh>

#include <aleph/utilities/Filesystem.hh>

#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <regex>
#include <vector>

#include <getopt.h>

#include <cmath>

using DataType           = double;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

struct Input
{
  std::string filename;
  PersistenceDiagram persistenceDiagram;

  std::string name;   // data set name
  unsigned dimension; // dimension of diagram
};

void usage()
{
  std::cerr << "Usage: persistence_diagram_statistics FILES\n"
            << "\n"
            << "Given a set of persistence diagrams, calculates numerous statistics\n"
            << "and writes them to STDOUT in CSV format.\n"
            << "\n"
            << "Currently, the following statistics are calculated:\n"
            << "  - Average persistence\n"
            << "  - Infinity norm\n"
            << "  - Total persistence\n"
            << "\n"
            << "Optional arguments:\n"
            << "\n"
            << " --invalid: Use the specified value to ignore certain values in every\n"
            << "            persistence diagram. This is useful if invalid values are\n"
            << "            encoded in the data.\n"
            << "\n"
            << " --power  : Use the specified power as an exponent during persistence\n"
            << "            calculations. This does not apply to the infinity norm of\n"
            << "            a persistence diagram.\n"
            << "\n\n";
}

std::pair<std::string, unsigned> parseFilename( std::string filename )
{
  if( aleph::utilities::extension( filename ) == ".txt" )
  {
    // Ignore anything else before the stem of the filename. Since we
    // are interested in parsing the dimension anyway, this is okay.
    filename = aleph::utilities::stem( filename );

    std::regex reDataSetPrefix( "(.*)_[dk]([[:digit:]]+)" );
    std::smatch matches;

    // Check if a recognizable prefix and suffix exist so that we may
    // grab information about the data set and its dimension. If not,
    // use the complete filename to identify the data set.
    if( std::regex_match( filename, matches, reDataSetPrefix ) )
    {
      std::string name   = matches[1];
      unsigned dimension = unsigned( std::stoul( matches[2] ) );

      return std::make_pair( name, dimension );
    }
  }

  // If we cannot parse the filename, just return its basename along
  // with dimension zero. This is a safe bet.
  return std::make_pair( aleph::utilities::stem( filename ), 0 );
}

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "invalid"       , required_argument, nullptr, 'i' },
    { "power"         , required_argument, nullptr, 'p' },
    { nullptr         , 0                , nullptr,  0  }
  };

  DataType invalid = std::numeric_limits<DataType>::has_quiet_NaN ? std::numeric_limits<DataType>::quiet_NaN() : std::numeric_limits<DataType>::max();
  double p         = 2.0;

  {
    int option = 0;
    while( ( option = getopt_long( argc, argv, "i:p:", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'i':
        invalid = static_cast<DataType>( std::stod( optarg ) );
        break;

      default:
        p = std::stod( optarg );
        break;
      }
    }
  }

  if( argc - optind < 1 )
  {
    usage();
    return -1;
  }

  std::vector<Input> inputs;
  inputs.reserve( static_cast<std::size_t>( argc - optind ) );

  std::vector<std::string> columns = {
    "file",
    "name",
    "dimension",
    "power",
    "total_persistence",
    "total_persistence_normalized",
    "infinity_norm",
    "average_persistence"
  };

  for( int i = optind; i < argc; i++ )
  {
    std::string filename = argv[i];

    std::cerr << "* Loading '" << filename << "'...";

    std::string name;
    unsigned dimension;

    std::tie( name, dimension ) = parseFilename( filename );

    Input input = {
      filename,
      aleph::io::load<DataType>( filename ),
      name,
      dimension
    };

    inputs.push_back( input );

    std::cerr << "finished\n";
  }

  // Header ------------------------------------------------------------

  {
    for( auto&& column : columns )
    {
      if( column != columns.front() )
        std::cout << ",";

      std::cout << column;
    }
    std::cout << "\n";
  }

  for( auto&& input : inputs )
  {
    if( !std::isnan( invalid ) && invalid != std::numeric_limits<DataType>::max() )
    {
      std::cerr << "* Filtering all persistence pairs that contain '" << invalid << "'...";

      using Point = typename PersistenceDiagram::Point;

      std::transform( input.persistenceDiagram.begin(), input.persistenceDiagram.end(),
                      input.persistenceDiagram.begin(),
                      [&invalid] ( const Point& p )
                      {
                        if( p.x() == invalid || p.y() == invalid )
                          return Point( DataType(), DataType() );
                        else
                          return Point( p );
                      } );

      input.persistenceDiagram.removeDiagonal();

      std::cerr << "\n";
    }

    auto totalPersistence           = aleph::totalPersistence( input.persistenceDiagram, p, false );
    auto totalPersistenceNormalized = totalPersistence / static_cast<decltype(totalPersistence)>( input.persistenceDiagram.size() );
    auto infinityNorm               = aleph::infinityNorm( input.persistenceDiagram );

    std::vector<DataType> persistence;
    aleph::persistence( input.persistenceDiagram, std::back_inserter( persistence ) );

    auto averagePersistence         = std::accumulate( persistence.begin(), persistence.end(), 0.0 ) / static_cast<double>( input.persistenceDiagram.size() );

    std::cout << "'" << input.filename      << "'" << ","
              << input.name                 << ","
              << input.dimension            << ","
              << p                          << ","
              << totalPersistence           << ","
              << totalPersistenceNormalized << ","
              << infinityNorm               << ","
              << averagePersistence         << "\n";
  }
}
