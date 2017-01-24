#include <iostream>
#include <string>
#include <vector>

#include "distances/Wasserstein.hh"

#include "persistenceDiagrams/IO.hh"
#include "persistenceDiagrams/PersistenceDiagram.hh"

using DataType           = double;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

int main( int argc, char** argv )
{
  if( argc <= 2 )
  {
    // TODO: Show usage
    return -1;
  }

  std::vector<std::string> filenames;
  filenames.reserve( argc - 2 );

  for( int i = 1; i < argc; i++ )
    filenames.push_back( argv[i] );

  std::vector<PersistenceDiagram> persistenceDiagrams;
  persistenceDiagrams.reserve( filenames.size() );

  for( auto&& filename : filenames )
  {
    std::cerr << "* Processing '" << filename << "'...";

    PersistenceDiagram persistenceDiagram = aleph::load<DataType>( filename );
    persistenceDiagrams.push_back( persistenceDiagram );

    std::cerr << "finished\n";
  }

  // TODO: Make configurable
  DataType power = DataType(1);

  for( std::size_t i = 0; i < persistenceDiagrams.size(); i++ )
  {
    auto&& pd1 = persistenceDiagrams.at(i);

    for( std::size_t j = i+1; j < persistenceDiagrams.size(); j++ )
    {
      auto&& pd2    = persistenceDiagrams.at(j);
      auto distance = aleph::distances::wassersteinDistance( pd1, pd2, power );

      std::cerr << "M[" << i << "," << j << "] = " << distance << "\n";
    }
  }
}
