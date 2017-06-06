/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  It analysis various aspects, such as the total persistence, of a set
  of persistence diagrams and writes all statistics to STDOUT.
*/

#include "persistenceDiagrams/Extraction.hh"
#include "persistenceDiagrams/Norms.hh"
#include "persistenceDiagrams/PersistenceDiagram.hh"

#include "persistenceDiagrams/io/Raw.hh"

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

void usage()
{
  // TODO: not yet implemented
}

int main( int argc, char** argv )
{
  if( argc <= 1 )
  {
    usage();
    return -1;
  }

  std::vector<Input> inputs;
  inputs.reserve( argc - 1 );

  std::vector<std::string> columns = {
    "file" ,
    "power",
    "total_persistence",
    "total_persistence_weighted",
    "infinity_norm"
  };

  for( int i = 1; i < argc; i++ )
  {
    std::string filename = argv[i];

    std::cerr << "* Loading '" << filename << "'...";

    Input input = {
      filename,
      aleph::io::load<DataType>( filename )
    };

    inputs.push_back( input );

    std::cerr << "finished\n";
  }

  // TODO: Make configurable
  double p = 2.0;

  for( auto&& column : columns )
    std::cout << column << " ";

  for( auto&& input : inputs )
  {
    auto totalPersistence         = aleph::totalPersistence( input.persistenceDiagram, p, false );
    auto totalPersistenceWeighted = aleph::totalPersistence( input.persistenceDiagram, p, true  );
    auto infinityNorm             = aleph::infinityNorm( input.persistenceDiagram );

    std::cout << "'" << input.filename    << "'"
              << " "
              << p
              << " "
              << totalPersistence         << " "
              << totalPersistenceWeighted << " "
              << infinityNorm             << "\n";
  }
}
