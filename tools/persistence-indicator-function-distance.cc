#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <string>
#include <vector>

#include <cmath>

// TODO: Replace this as soon as possible with a more modern option
// parser interface.
#include <getopt.h>

#include "distances/Wasserstein.hh"

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
  std::string name;
  std::string filename;
  unsigned dimension;

  PersistenceDiagram persistenceDiagram;
  PersistenceIndicatorFunction persistenceIndicatorFunction;
};

void storeMatrix( const std::vector< std::vector<double> >& M, std::ostream& out )
{
  if( M.empty() )
    return;

  auto rows = M.size();
  auto cols = M.front().size();

  for( decltype(rows) row = 0; row < rows; row++ )
  {
    for( decltype(cols) col = 0; col < cols; col++ )
    {
      if( col != 0 )
        out << " ";

      out << M[row][col];
    }

    out << "\n";
  }
}

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
      return PersistenceIndicatorFunction( it->persistenceIndicatorFunction );
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

double wassersteinDistance( const std::vector<DataSet>& dataSet1, const std::vector<DataSet>& dataSet2, unsigned minDimension, unsigned maxDimension, double power )
{
  auto getPersistenceDiagram = [] ( const std::vector<DataSet>& dataSet, unsigned dimension )
  {
    auto it = std::find_if( dataSet.begin(), dataSet.end(),
                            [&dimension] ( const DataSet& dataSet )
                            {
                              return dataSet.dimension == dimension;
                            } );

    if( it != dataSet.end() )
      return PersistenceDiagram( it->persistenceDiagram );
    else
      return PersistenceDiagram();
  };

  double d = 0.0;

  for( unsigned dimension = minDimension; dimension <= maxDimension; dimension++ )
  {
    auto D1 = getPersistenceDiagram( dataSet1, dimension );
    auto D2 = getPersistenceDiagram( dataSet2, dimension );

    d += aleph::distances::wassersteinDistance( D1, D2, power );
  }

  d = std::pow( d, 1.0 / power );
  return d;
}


int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "power"      , required_argument, nullptr, 'p' },
    { "wasserstein", no_argument      , nullptr, 'w' },
    { nullptr      , 0                , nullptr,  0  }
  };

  double power                = 2.0;
  bool useWassersteinDistance = false;

  int option = 0;
  while( ( option = getopt_long( argc, argv, "p:w", commandLineOptions, nullptr ) ) != -1 )
  {
    switch( option )
    {
    case 'p':
      power = std::stod( optarg );
      break;
    case 'w':
      useWassersteinDistance = true;
      break;
    default:
      break;
    }
  }

  if( ( argc - optind ) <= 1 )
  {
    // TODO: Show usage
    return -1;
  }

  // Maps filenames to indices. I need this to ensure that the internal
  // ordering of files coincides with the shell's ordering.
  std::map<std::string, unsigned> filenameMap;

  std::vector< std::vector<DataSet> > dataSets;

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

    unsigned index = 0;

    for( int i = optind; i < argc; i++ )
    {
      filenames.push_back( argv[i] );

      if( std::regex_match( filenames.back(), matches, reDataSetPrefix ) )
      {
        auto name = matches[1];

        if( filenameMap.find( name ) == filenameMap.end() )
          filenameMap[ name ] = index++;
      }
    }

    dataSets.resize( filenameMap.size() );

    for( auto&& filename : filenames )
    {
      if( std::regex_match( filename, matches, reDataSetPrefix ) )
      {
        auto name      = matches[1];
        auto dimension = unsigned( std::stoul( matches[2] ) );

        dataSets.at( filenameMap[name] ).push_back( { name, filename, dimension, {}, {} } );

        minDimension = std::min( minDimension, dimension );
        maxDimension = std::max( maxDimension, dimension );
      }
    }
  }

  // Load persistence diagrams & calculate indicator functions ---------

  std::vector<PersistenceDiagram> persistenceDiagrams;

  for( auto&& sets : dataSets )
  {
    for( auto&& dataSet : sets )
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
    for( auto it2 = dataSets.begin(); it2 != dataSets.end(); ++it2, ++col )
    {
      if( it1 == it2 )
      {
        distances[row][row] = 0.0;
        continue;
      }

      double d = 0.0;

      if( useWassersteinDistance )
        d = wassersteinDistance( *it1, *it2, minDimension, maxDimension, power );
      else
        d = distance( *it1, *it2, minDimension, maxDimension );

      distances[row][col] = d;
      distances[col][row] = d;
    }
  }

  std::cerr << "Storing matrix...";

  storeMatrix( distances, std::cout );

  std::cerr << "finished\n";
}
