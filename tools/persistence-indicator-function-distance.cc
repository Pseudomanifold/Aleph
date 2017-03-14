#include <iostream>
#include <map>
#include <regex>
#include <string>
#include <vector>

#include "persistenceDiagrams/IO.hh"
#include "persistenceDiagrams/PersistenceDiagram.hh"
#include "persistenceDiagrams/PersistenceIndicatorFunction.hh"

#include "utilities/Filesystem.hh"

using DataType                     = double;
using PersistenceDiagram           = aleph::PersistenceDiagram<DataType>;
using PersistenceIndicatorFunction = aleph::math::StepFunction<DataType>;

/*
  Auxiliary structure for describing a data set. I need this in order to
  figure out the corresponding dimension of the persistence diagram.
*/

struct DataSet
{
  std::string filename;
  unsigned dimension;

  PersistenceDiagram persistenceDiagram;
  PersistenceIndicatorFunction persistenceIndicatorFunction;
};

int main( int argc, char** argv )
{
  if( argc <= 1 )
  {
    // TODO: Show usage
    return -1;
  }

  // Maps data sets to file names. Whether a new filename constitutes
  // a new data set is decided by taking a look at the prefix of the
  // filename.
  std::map<std::string, std::vector<DataSet> > dataSets;

  // TODO: Make this regular expression more, well, 'expressive' and
  // support more methods of specifying a dimension.
  std::regex reDataSetPrefix( "(.*)_k([[:digit:]]+)\\.txt" );
  std::smatch matches;

  // Get filenames & prefixes ------------------------------------------

  {
    std::vector<std::string> filenames;
    filenames.reserve( argc - 1 );

    for( int i = 1; i < argc; i++ )
      filenames.push_back( argv[i] );

    for( auto&& filename : filenames )
    {
      if( std::regex_match( filename, matches, reDataSetPrefix ) )
      {
        auto name      = matches[1];
        auto dimension = unsigned( std::stoul( matches[2] ) );

        dataSets[name].push_back( { filename, dimension, {}, {} } );
      }
    }
  }

  // Load persistence diagrams -----------------------------------------

  std::vector<PersistenceDiagram> persistenceDiagrams;

  for( auto&& pair : dataSets )
  {
    for( auto&& dataSet : pair.second )
    {
      std::cerr << "* Processing '" << dataSet.filename << "'...";

      dataSet.persistenceDiagram = aleph::load<DataType>( dataSet.filename );

      // FIXME: This is only required in order to ensure that the
      // persistence indicator function has a finite integral; it
      // can be solved more elegantly by using a special value to
      // indicate infinite intervals.
      dataSet.persistenceDiagram.removeUnpaired();

      dataSet.persistenceIndicatorFunction
         = aleph::persistenceIndicatorFunction( dataSet.persistenceDiagram );

      std::cerr << "finished\n";
    }
  }
}
