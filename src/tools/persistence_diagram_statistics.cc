/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  It analyses various aspects, such as the total persistence, of a set
  of persistence diagrams and writes all statistics to STDOUT.
*/

#include <aleph/persistenceDiagrams/Extraction.hh>
#include <aleph/persistenceDiagrams/Norms.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistenceDiagrams/io/Raw.hh>

#include <iostream>
#include <numeric>
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
    "infinity_norm",
    "average_persistence",
    "average_persistence_weighted"
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

  // TODO: Add explanatory header with comment symbols?
  for( auto&& column : columns )
    std::cout << column << " ";
  std::cout << "\n";

  for( auto&& input : inputs )
  {
    auto totalPersistence         = aleph::totalPersistence( input.persistenceDiagram, p, false );
    auto totalPersistenceWeighted = aleph::totalPersistence( input.persistenceDiagram, p, true  );
    auto infinityNorm             = aleph::infinityNorm( input.persistenceDiagram );

    std::vector<DataType> persistence;
    aleph::persistence( input.persistenceDiagram, std::back_inserter( persistence ) );

    std::vector<DataType> weightedPersistence;
    aleph::weightedPersistence( input.persistenceDiagram, std::back_inserter( weightedPersistence ) );

    auto averagePersistence         = std::accumulate( persistence.begin(), persistence.end(), 0.0 ) / static_cast<double>( input.persistenceDiagram.size() );
    auto averageWeightedPersistence = std::accumulate( weightedPersistence.begin(), weightedPersistence.end(), 0.0 ) / static_cast<double>( input.persistenceDiagram.size() );

    std::cout << "'" << input.filename      << "'" << " "
              << p                          << " "
              << totalPersistence           << " "
              << totalPersistenceWeighted   << " "
              << infinityNorm               << " "
              << averagePersistence         << " "
              << averageWeightedPersistence << "\n";
  }
}
