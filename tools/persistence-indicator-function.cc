#include <iostream>
#include <string>
#include <vector>

#include "persistenceDiagrams/IO.hh"
#include "persistenceDiagrams/PersistenceDiagram.hh"
#include "persistenceDiagrams/PersistenceIndicatorFunction.hh"

using DataType           = double;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

int main( int argc, char** argv )
{
  if( argc <= 2 )
  {
    // TODO: Show usage
    return -1;
  }

  // Get filenames -----------------------------------------------------

  std::vector<std::string> filenames;
  filenames.reserve( argc - 2 );

  for( int i = 1; i < argc; i++ )
    filenames.push_back( argv[i] );

  // Load persistence diagrams -----------------------------------------

  std::vector<PersistenceDiagram> persistenceDiagrams;
  persistenceDiagrams.reserve( filenames.size() );

  for( auto&& filename : filenames )
  {
    std::cerr << "* Processing '" << filename << "'...";

    PersistenceDiagram persistenceDiagram = aleph::load<DataType>( filename );
    persistenceDiagrams.push_back( persistenceDiagram );

    std::cerr << "finished\n";
  }

  // Calculate persistence indicator functions -------------------------

  for( auto&& D : persistenceDiagrams )
    aleph::persistenceIndicatorFunction( D );
}
