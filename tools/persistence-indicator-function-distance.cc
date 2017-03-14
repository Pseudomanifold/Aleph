#include <algorithm>
#include <iostream>
#include <limits>
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

// Calculates the distance between two data sets. This requires enumerating all
// dimensions, and finding the corresponding persistence indicator function. If
// no such function exists, the function uses the norm of the function.
double distance( const std::vector<DataSet>& dataSet1, const std::vector<DataSet>& dataSet2, unsigned minDimension, unsigned maxDimension )
{
  auto getPersistenceIndicatorFunction = [] ( const std::vector<DataSet>& dataSet, unsigned dimension )
  {
    auto it = std::find_if( dataSet.begin(), dataSet.end(),
                            [&dimension] ( const DataSet& dataSet )
                            {
                              return dataSet.dimension == dimension;
                            } );

    if( it != dataSet.end() )
      return it->persistenceIndicatorFunction;
    else
      return PersistenceIndicatorFunction();
  };

  double d = 0.0;

  for( unsigned dimension = minDimension; dimension <= maxDimension; dimension++ )
  {
    auto f = getPersistenceIndicatorFunction( dataSet1, dimension );
    auto g = getPersistenceIndicatorFunction( dataSet2, dimension );

    g = g * (-1.0);
    d = d + (f+g).integral();
  }

  return d;
}

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

  unsigned minDimension = std::numeric_limits<unsigned>::max();
  unsigned maxDimension = 0;

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

        minDimension = std::min( minDimension, dimension );
        maxDimension = std::max( maxDimension, dimension );
      }
    }
  }

  // Load persistence diagrams & calculate indicator functions ---------

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

  // Calculate all distances -------------------------------------------

  std::vector< std::vector<double> > distances;
  distances.resize( dataSets.size(), std::vector<double>( dataSets.size() ) );

  std::size_t row = 0;

  for( auto it1 = dataSets.begin(); it1 != dataSets.end(); ++it1, ++row )
  {
    std::size_t col = 0;
    for( auto it2 = std::next( it1 ); it2 != dataSets.end(); ++it2, ++col )
    {
      auto d = distance( it1->second, it2->second, minDimension, maxDimension );

      distances[row][col] = d;
      distances[col][row] = d;

      std::cerr << "D[" << row << "][" << col << "]: " << d << "\n";
    }
  }
}
