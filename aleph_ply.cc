#include "io/PLY.hh"

#include <iostream>

using DataType   = double;
using VertexType = unsigned;

int main( int argc, char** argv )
{
  std::string filename;
  std::string property = "quality";

  if( argc == 1 )
    return -1;

  if( argc >= 2 )
    filename = argv[1];

  if( argc >= 3 )
    property = argv[2];

  auto K = aleph::io::loadPLY<DataType, VertexType>( filename, property );

  std::cout << "* Loaded simplicial complex with " << K.size() << " simplices\n";
}
