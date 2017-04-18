#include "persistenceDiagrams/IO.hh"
#include "persistenceDiagrams/Norms.hh"
#include "persistenceDiagrams/PersistenceDiagram.hh"

#include <iostream>
#include <string>
#include <vector>

using DataType           = double;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

struct Input
{
  std::string filename;
  PersistenceDiagram persistenceDiagram;
};

int main( int argc, char** argv )
{
  std::vector<Input> inputs;
  inputs.reserve( argc - 1 );

  for( int i = 1; i < argc; i++ )
  {
    std::string filename = argv[i];

    std::cerr << "* Loading '" << filename << "'...";

    Input input = {
      filename,
      aleph::load<DataType>( filename )
    };

    inputs.push_back( input );

    std::cerr << "finished\n";
  }

  for( auto&& input : inputs )
  {
    // TODO: Make $p$ user-selectable
    std::cout << "# " << input.filename << "\n"
              << "Total persistence: " << aleph::totalPersistence( input.persistenceDiagram ) << "\n"
              << "p-norm:            " << aleph::pNorm( input.persistenceDiagram )            << "\n"
              << "p:                 " << 2.0                                                 << "\n";
  }
}
