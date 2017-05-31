#include <algorithm>
#include <functional>
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

#include "distances/Hausdorff.hh"
#include "distances/Wasserstein.hh"

#include "persistenceDiagrams/PersistenceDiagram.hh"
#include "persistenceDiagrams/PersistenceIndicatorFunction.hh"

#include "persistenceDiagrams/io/Raw.hh"

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

/*
  Stores a matrix in an output stream. The matrix is formatted such that
  individual values are separated by spaces and each row ends with '\n'.

  This format can be easily parsed by auxiliary programs such as gnuplot
  or R.
*/

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

/*
  Calculates the topological distance between two data sets using persistence
  indicator functions. This requires enumerating all dimensions and finding a
  corresponding persistence indicator function. If no suitable function could
  be found, the calculation defaults to calculating the norm.
*/

double distancePIF( const std::vector<DataSet>& dataSet1,
                    const std::vector<DataSet>& dataSet2,
                    unsigned minDimension,
                    unsigned maxDimension,
                    double power )
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

    g = -g;
    if( power == 1.0 )
      d = d + (f+g).abs().integral();
    else
      d = d + (f+g).abs().integral_p( power );
  }

  return d;
}

/*
  Calculates the topological distance between two data sets, using
  a standard distance between two persistence diagrams, for example
  the Hausdorff, Wasserstein, or bottleneck distance.

  By default, the Wasserstein distance is calculated.
*/

template <class Functor>
double persistenceDiagramDistance( const std::vector<DataSet>& dataSet1,
                                   const std::vector<DataSet>& dataSet2,
                                   unsigned minDimension,
                                   unsigned maxDimension,
                                   double power,
                                   Functor functor = [] ( const PersistenceDiagram& D1, const PersistenceDiagram& D2, double power )
                                   {
                                     return aleph::distances::wassersteinDistance( D1, D2, power );
                                   } )
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

    d += functor( D1, D2, power );
  }

  d = std::pow( d, 1.0 / power );
  return d;
}


int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "power"      , required_argument, nullptr, 'p' },
    { "exp"        , no_argument      , nullptr, 'e' },
    { "hausdorff"  , no_argument      , nullptr, 'h' },
    { "indicator"  , no_argument      , nullptr, 'i' },
    { "kernel"     , no_argument      , nullptr, 'k' },
    { "wasserstein", no_argument      , nullptr, 'w' },
    { nullptr      , 0                , nullptr,  0  }
  };

  double power                      = 2.0;
  bool useExponentialFunction       = false;
  bool useIndicatorFunctionDistance = false;
  bool calculateKernel              = false;
  bool useWassersteinDistance       = false;

  int option = 0;
  while( ( option = getopt_long( argc, argv, "p:ehikw", commandLineOptions, nullptr ) ) != -1 )
  {
    switch( option )
    {
    case 'p':
      power = std::stod( optarg );
      break;
    case 'e':
      useExponentialFunction = true;
      break;
    case 'h':
      useWassersteinDistance       = false;
      useIndicatorFunctionDistance = false;
      break;
    case 'i':
      useIndicatorFunctionDistance = true;
      useWassersteinDistance       = false;
      break;
    case 'k':
      calculateKernel = true;
      break;
    case 'w':
      useIndicatorFunctionDistance = false;
      useWassersteinDistance       = true;
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
  std::regex reDataSetPrefix( "(.*)_[dk]([[:digit:]]+)\\.txt" );
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

      auto name = filenames.back();

      // Check whether the name contains a recognizable prefix and
      // suffix. If not, use the complete filename to identify the
      // data set.
      if( std::regex_match( filenames.back(), matches, reDataSetPrefix ) )
        name = matches[1];

      if( filenameMap.find( name ) == filenameMap.end() )
        filenameMap[ name ] = index++;
    }

    dataSets.resize( filenameMap.size() );

    for( auto&& filename : filenames )
    {
      auto name      = filename;
      auto dimension = 0u;

      // Check if a recognizable prefix and suffix exist so that we may
      // grab information about the data set and its dimension. If not,
      // use the complete filename to identify the data set.
      if( std::regex_match( filename, matches, reDataSetPrefix ) )
      {
        name      = matches[1];
        dimension = unsigned( std::stoul( matches[2] ) );
      }

      dataSets.at( filenameMap[name] ).push_back( { name, filename, dimension, {}, {} } );

      minDimension = std::min( minDimension, dimension );
      maxDimension = std::max( maxDimension, dimension );
    }
  }

  // Load persistence diagrams & calculate indicator functions ---------

  std::vector<PersistenceDiagram> persistenceDiagrams;

  for( auto&& sets : dataSets )
  {
    for( auto&& dataSet : sets )
    {
      std::cerr << "* Processing '" << dataSet.filename << "'...";

      dataSet.persistenceDiagram = aleph::io::load<DataType>( dataSet.filename );

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

  // Setup distance functor --------------------------------------------

  std::function< double( const PersistenceDiagram&, const PersistenceDiagram&, double ) > functor
    = !useIndicatorFunctionDistance ? useWassersteinDistance ? [] ( const PersistenceDiagram& D1, const PersistenceDiagram& D2, double p )
                                                               {
                                                                 return aleph::distances::wassersteinDistance( D1, D2, p );
                                                               }
                                                             : [] ( const PersistenceDiagram& D1, const PersistenceDiagram& D2, double p )
                                                               {
                                                                 return std::pow( aleph::distances::hausdorffDistance( D1, D2 ), p );
                                                               }
                                    : [] ( const PersistenceDiagram&, const PersistenceDiagram&, double )
                                    {
                                      return 0.0;
                                    };

  // Calculate all distances -------------------------------------------

  std::vector< std::vector<double> > distances;
  distances.resize( dataSets.size(), std::vector<double>( dataSets.size() ) );

  #pragma omp parallel for collapse(2)
  for( std::size_t row = 0; row < dataSets.size(); row++ )
  {
    for( std::size_t col = 0; col < dataSets.size(); col++ )
    {
      if( row <= col )
        continue;

      double d = 0.0;

      if( useIndicatorFunctionDistance )
        d = distancePIF( dataSets.at(row), dataSets.at(col), minDimension, maxDimension, power );
      else
        d = persistenceDiagramDistance( dataSets.at(row), dataSets.at(col), minDimension, maxDimension, power, functor );

      if( calculateKernel )
      {
        d = -d;
        if( useExponentialFunction )
          d = std::exp( d );
      }

      distances[row][col] = d;
      distances[col][row] = d;
    }
  }

  std::cerr << "Storing matrix...";

  storeMatrix( distances, std::cout );

  std::cerr << "finished\n";
}
