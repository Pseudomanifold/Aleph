#include <aleph/topology/filtrations/Data.hh>

#include <aleph/topology/io/Function.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/utilities/String.hh>

#include <fstream>
#include <istream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <getopt.h>

using DataType          = double;
using VertexType        = unsigned;
using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

std::vector<SimplicialComplex> readData( std::istream& in, bool useSublevelSetFiltration )
{
  std::vector<SimplicialComplex> complexes;
  std::string line;

  while( std::getline( in, line ) )
  {
    auto tokens
      = aleph::utilities::split(
          line,
          std::string( "[:;,[:space:]]+" )
    );

    std::vector<DataType> values;
    values.reserve( tokens.size() );

    for( auto&& token : tokens )
    {
      bool success = false;
      auto value   = aleph::utilities::convert<DataType>( token, success );

      if( !success )
        throw std::runtime_error( "Unable to convert token to expected data type" );

      values.emplace_back( value );
    }

    complexes.push_back(
      aleph::topology::io::loadFunction<SimplicialComplex>(
        values.begin(), values.end(),
        [&useSublevelSetFiltration] ( DataType x, DataType y )
        {
          if( useSublevelSetFiltration )
            return std::max(x,y);
          else
            return std::min(x,y);
        }
      )
    );
  }

  return complexes;
}

void usage()
{
}

int main( int argc, char** argv )
{
  // Options parsing ---------------------------------------------------
  //
  // By default, a sublevel set filtration is being calculated for the
  // input data set. One or more input data sets may be specified at a
  // time. Using '-' indicates that input should be read from `stdin`.

  bool useSublevelSetFiltration = true;

  {
    static option commandLineOptions[] =
    {
      { "sublevels"  , no_argument, nullptr, 's' },
      { "superlevels", no_argument, nullptr, 'S' },
      { nullptr      , 0          , nullptr,  0  }
    };

    int option = 0;
    while( ( option = getopt_long( argc, argv, "sS", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 's':
        useSublevelSetFiltration = true;
        break;

      case 'S':
        useSublevelSetFiltration = false;
        break;

      default:
        break;
      }
    }
  }

  if( ( argc - optind ) < 1 )
  {
    usage();
    return -1;
  }

  std::vector<SimplicialComplex> complexes;

  for( int i = optind; i < argc; i++ )
  {
    std::string filename = argv[i];
    std::cerr << "* Reading '" << filename << "'...";

    std::istream* in = &std::cin;
    std::ifstream fin;

    if( filename != "-" and not filename.empty() )
    {
      fin.open( filename );
      in = &fin;
    }

    auto localComplexes
      = readData( *in, useSublevelSetFiltration );

    std::cerr << "finished\n";

    complexes.insert( complexes.end(),
      localComplexes.begin(), localComplexes.end()
    );
  }

  std::cerr << "* Read " << complexes.size() << " simplicial complexes\n";
}
