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

using DataType          = double;
using VertexType        = unsigned;
using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

std::vector<SimplicialComplex> readData( std::istream& in )
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
        [] ( DataType x, DataType y )
        {
          return std::max(x,y);
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
  if( argc < 1 )
  {
    usage();
    return -1;
  }
}
