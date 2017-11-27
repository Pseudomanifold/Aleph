#include <aleph/math/Bootstrap.hh>
#include <aleph/math/PiecewiseLinearFunction.hh>

#include <aleph/utilities/String.hh>

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

#include <getopt.h>

using DataType = double;
using Function = aleph::math::PiecewiseLinearFunction<DataType>;

Function load( const std::string& filename )
{
  std::ifstream in( filename );
  if( !in )
    throw std::runtime_error( "Unable to open file for reading" );

  std::string line;
  std::vector< std::pair<DataType, DataType> > data;

  while( std::getline( in, line ) )
  {
    line = aleph::utilities::trim( line );

    if( line.empty() || line.front() == '#' )
      continue;

    DataType x = DataType();
    DataType y = DataType();

    std::istringstream converter( line );

    converter >> x >> y;

    if( !converter )
      throw std::runtime_error( "Conversion failed" );

    data.push_back( std::make_pair(x,y) );
  }

  return Function( data.begin(), data.end() );
}

int main( int argc, char** argv )
{
  DataType alpha               = 0.95;
  unsigned numBootstrapSamples = 10;

  {
    static option commandLineOptions[] =
    {
      { "alpha", required_argument, nullptr, 'a' },
      { "n"    , required_argument, nullptr, 'n' },
      { nullptr, 0                , nullptr,  0  }
    };

    int option = 0;
    while( ( option = getopt_long( argc, argv, "a:", commandLineOptions, nullptr ) ) != - 1 )
    {
      switch( option )
      {
      case 'a':
        alpha = DataType( std::stod( optarg ) );
        break;
      case 'n':
        numBootstrapSamples = unsigned( std::stoull( optarg ) );
        break;
      default:
        break;
      }
    }
  }

  if( argc - optind <= 1 )
    return -1;

  // Load functions ----------------------------------------------------

  std::vector<Function> functions;

  for( int i = optind; i < argc; i++ )
    functions.push_back( load( argv[i] ) );
}
