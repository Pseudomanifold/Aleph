/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  Given a set of persistence diagrams, it calculates various topological
  distances and returns a distance matrix.

  This tool can be helpful in different application scenarios:

  1. You want to determine the dissimilarity between two
     high-dimensional shapes, filtered by their distance
     function.

  2. You want to measure how a data descriptor, e.g. any
     density estimator, is changing over embeddings of a
     high-dimensional data set.

  3. You want to determine if certain samples of a space
     have the same characteristics than the original.

  The tool attempts to be smart and groups different inputs according to
  their common prefix. Currently, it only understands _d and _k as valid
  suffixes. Hence, the following input files are considered to belong to
  the same data set:

  - Test_d01
  - Test_d05
  - Test_d07

  Likewise:

  - Test_k1
  - Test_k7
  - Test_k9

  Please keep this in mind when using the tool.
*/

#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <string>
#include <vector>

#include <cmath>

#include <getopt.h>

#include <aleph/persistenceDiagrams/Envelope.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/PersistenceIndicatorFunction.hh>

#include <aleph/persistenceDiagrams/distances/Hausdorff.hh>
#include <aleph/persistenceDiagrams/distances/Wasserstein.hh>

#include <aleph/persistenceDiagrams/io/JSON.hh>
#include <aleph/persistenceDiagrams/io/Raw.hh>

#include <aleph/utilities/Filesystem.hh>

using DataType                     = double;
using PersistenceDiagram           = aleph::PersistenceDiagram<DataType>;
using PersistenceIndicatorFunction = aleph::math::StepFunction<DataType>;
using EnvelopeFunction             = aleph::math::PiecewiseLinearFunction<DataType>;

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
  EnvelopeFunction envelopeFunction;
};

/* Usage information */
void usage()
{
  std::cerr << "Usage: topological_distance [--power=POWER] [--kernel] [--exp] [--sigma]\n"
            << "                            [--hausdorff|envelope|indicator|wasserstein]\n"
            << "                            [--clean] [--factor=FACTOR] FILES\n"
            << "\n"
            << "Calculates distances between a set of persistence diagrams, stored\n"
            << "in FILES. By default, this tool calculates Hausdorff distances for\n"
            << "all diagrams. This can be modified.\n"
            << "\n"
            << "If no other value is given, all distances are weighted using $p=2$\n"
            << "during the construction of a pairwise distance matrix. Furthermore\n"
            << "this tool can calculate kernels for use in kernel-based methods in\n"
            << "machine learning.\n"
            << "\n"
            << "Use --factor=FACTOR to specify the factor that will be used in the\n"
            << "treatment of unpaired points. If set to any non-zero value, all of\n"
            << "the unpaired points are multiplied by it.\n"
            << "\n"
            << "The distance matrix is written to STDOUT. Rows and columns will be\n"
            << "separated by whitespace.\n"
            << "\n"
            << "This tool tries to be smart and is able to detect whether a set of\n"
            << "persistence diagrams belongs to the same group. This works only if\n"
            << "each file contains a suffix with digits that is preceded by either\n"
            << "a 'd' (for dimension) or a 'k' (for clique dimension).\n"
            << "\n"
            << "Flags:\n"
            << "  -c: clean persistence diagrams (remove unpaired points)\n"
            << "  -e: use exponential weighting for kernel calculation\n"
            << "  -E: calculate envelope function distances\n"
            << "  -h: calculate Hausdorff distances\n"
            << "  -i: calculate persistence indicator function distances\n"
            << "  -k: calculate kernel values instead of distances\n"
            << "  -n: normalize the persistence indicator function\n"
            << "  -s: use sigma as a scale parameter for the kernel\n"
            << "  -w: calculate Wasserstein distances\n"
            << "\n";
}

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
                    double power,
                    bool normalize )
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

    if( normalize )
    {
      f = aleph::math::normalize( f );
      g = aleph::math::normalize( g );
    }

    g = -g;
    if( power == 1.0 )
      d = d + (f+g).abs().integral();
    else
      d = d + (f+g).abs().pow( power ).integral();
  }

  return d;
}

/*
  Calculates the topological distance between two data sets using
  envelope functions. This requires enumerating all dimensions in
  order to find the corresponding envelope function. If no proper
  function is found, the method defaults to calculating the norm.
*/

double distanceEnvelopeFunctions( const std::vector<DataSet>& dataSet1,
                                  const std::vector<DataSet>& dataSet2,
                                  unsigned minDimension,
                                  unsigned maxDimension,
                                  double power )
{
  auto getEnvelopeFunction = [] ( const std::vector<DataSet>& dataSet, unsigned dimension )
  {
    auto it = std::find_if( dataSet.begin(), dataSet.end(),
                            [&dimension] ( const DataSet& dataSet )
                            {
                              return dataSet.dimension == dimension;
                            } );

    if( it != dataSet.end() )
      return EnvelopeFunction( it->envelopeFunction );
    else
      return EnvelopeFunction();
  };

  double d = 0.0;

  for( unsigned dimension = minDimension; dimension <= maxDimension; dimension++ )
  {
    auto f = getEnvelopeFunction( dataSet1, dimension );
    auto g = getEnvelopeFunction( dataSet2, dimension );

    g = -g;
    if( power == 1.0 )
      d = d + (f+g).abs().integral();
    else
      d = d + (f+g).abs().integral( power );
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

DataType getMaximum( const PersistenceDiagram& diagram )
{
  DataType max = std::numeric_limits<DataType>::lowest();

  for( auto&& p : diagram )
  {
    max = std::max( max, p.x() );
    if( !p.isUnpaired() )
      max = std::max( max, p.y() );
  }

  return max;
}

PersistenceDiagram postprocess( const PersistenceDiagram& diagram, bool clean, DataType infinityFactor )
{
  auto result = diagram;

  if( clean )
  {
    result.removeDiagonal();
    result.removeUnpaired();
  }

  if( infinityFactor != DataType() )
  {
    auto max = getMaximum( result );

    std::transform( result.begin(), result.end(), result.begin(),
                    [&infinityFactor, &max] ( const PersistenceDiagram::Point& p )
                    {
                      if( p.isUnpaired() )
                        return PersistenceDiagram::Point( p.x(), infinityFactor * max );
                      else
                        return PersistenceDiagram::Point( p );
                    } );
  }

  return result;
}

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "factor"     , required_argument, nullptr, 'f' },
    { "power"      , required_argument, nullptr, 'p' },
    { "sigma"      , required_argument, nullptr, 's' },
    { "clean"      , no_argument      , nullptr, 'c' },
    { "envelope"   , no_argument      , nullptr, 'E' },
    { "exp"        , no_argument      , nullptr, 'e' },
    { "hausdorff"  , no_argument      , nullptr, 'h' },
    { "indicator"  , no_argument      , nullptr, 'i' },
    { "normalize"  , no_argument      , nullptr, 'n' },
    { "kernel"     , no_argument      , nullptr, 'k' },
    { "wasserstein", no_argument      , nullptr, 'w' },
    { nullptr      , 0                , nullptr,  0  }
  };

  DataType infinityFactor           = DataType();
  double power                      = 2.0;
  double sigma                      = 1.0;
  bool cleanPersistenceDiagrams     = false;
  bool useExponentialFunction       = false;
  bool useEnvelopeFunctionDistance  = false;
  bool useIndicatorFunctionDistance = false;
  bool normalize                    = false;
  bool calculateKernel              = false;
  bool useWassersteinDistance       = false;

  int option = 0;
  while( ( option = getopt_long( argc, argv, "f:p:s:ceEhinkw", commandLineOptions, nullptr ) ) != -1 )
  {
    switch( option )
    {
    case 'f':
      infinityFactor = static_cast<DataType>( std::stod( optarg ) );
      break;
    case 'p':
      power = std::stod( optarg );
      break;
    case 's':
      sigma = std::stod( optarg );
      break;
    case 'c':
      cleanPersistenceDiagrams = true;
      break;
    case 'E':
      useEnvelopeFunctionDistance  = true;
      useWassersteinDistance       = false;
      useIndicatorFunctionDistance = false;
      break;
    case 'e':
      useExponentialFunction = true;
      break;
    case 'h':
      useWassersteinDistance       = false;
      useIndicatorFunctionDistance = false;
      useEnvelopeFunctionDistance  = false;
      break;
    case 'i':
      useIndicatorFunctionDistance = true;
      useEnvelopeFunctionDistance  = false;
      useWassersteinDistance       = false;
      break;
    case 'k':
      calculateKernel = true;
      break;
    case 'n':
      normalize = true;
      break;
    case 'w':
      useEnvelopeFunctionDistance  = false;
      useIndicatorFunctionDistance = false;
      useWassersteinDistance       = true;
      break;
    default:
      break;
    }
  }

  if( ( argc - optind ) <= 1 )
  {
    usage();
    return -1;
  }

  std::vector< std::vector<DataSet> > dataSets;

  // Get filenames & prefixes ------------------------------------------

  unsigned minDimension = std::numeric_limits<unsigned>::max();
  unsigned maxDimension = 0;

  {
    std::vector<std::string> filenames;
    filenames.reserve( argc - 1 );

    for( int i = optind; i < argc; i++ )
      filenames.push_back( argv[i] );

    // This should never happen...
    if( filenames.empty() )
      return -1;

    // If the first filename is a text file, I am assuming that the rest
    // of them also are. The program will then read all diagrams, try to
    // match them to a dimension, and store them.
    if( aleph::utilities::extension( filenames.front() ) == ".txt" )
    {
      std::regex reDataSetPrefix( "(.*)_[dk]([[:digit:]]+)\\.txt" );
      std::smatch matches;

      unsigned index = 0;

      // Maps filenames to indices. I need this to ensure that the internal
      // ordering of files coincides with the shell's ordering.
      std::map<std::string, unsigned> filenameMap;

      for( auto&& filename : filenames )
      {
        auto name = filename;

        // Check whether the name contains a recognizable prefix and
        // suffix. If not, use the complete filename to identify the
        // data set.
        if( std::regex_match( filename, matches, reDataSetPrefix ) )
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

        dataSets.at( filenameMap[name] ).push_back( { name, filename, dimension, {}, {}, {} } );

        minDimension = std::min( minDimension, dimension );
        maxDimension = std::max( maxDimension, dimension );
      }

      // Load persistence diagrams & calculate indicator functions -----

      for( auto&& sets : dataSets )
      {
        for( auto&& dataSet : sets )
        {
          std::cerr << "* Processing '" << dataSet.filename << "'...";

          dataSet.persistenceDiagram
            = postprocess( aleph::io::load<DataType>( dataSet.filename ),
                           cleanPersistenceDiagrams,
                           infinityFactor );

          // FIXME: This is only required in order to ensure that the
          // persistence indicator function has a finite integral; it
          // can be solved more elegantly by using a special value to
          // indicate infinite intervals.
          auto pd = dataSet.persistenceDiagram;
          pd.removeUnpaired();

          dataSet.persistenceIndicatorFunction = aleph::persistenceIndicatorFunction( pd );
          dataSet.envelopeFunction             = aleph::Envelope()( pd );

          std::cerr << "finished\n";
        }
      }
    }
    else if( aleph::utilities::extension( filenames.front() ) == ".json" )
    {
      dataSets.reserve( filenames.size() );

      for( auto&& filename : filenames )
      {
        auto persistenceDiagrams= aleph::io::readJSON<DataType>( filename );

        std::vector<DataSet> dataSet;
        dataSet.reserve( persistenceDiagrams.size() );

        for( auto&& diagram : persistenceDiagrams )
        {
          diagram
            = postprocess( diagram,
                           cleanPersistenceDiagrams,
                           infinityFactor );

          auto dimension = static_cast<unsigned>( diagram.dimension() );
          minDimension   = std::min( minDimension, dimension );
          maxDimension   = std::max( maxDimension, dimension );

          auto name  = aleph::utilities::stem( filename );
          name      += "_";
          name      += "d" + std::to_string( diagram.dimension() );

          // FIXME: This is only required in order to ensure that the
          // persistence indicator function has a finite integral; it
          // can be solved more elegantly by using a special value to
          // indicate infinite intervals.
          auto pd = diagram;
          pd.removeUnpaired();

          dataSet.push_back( { name,
                               filename,
                               dimension,
                               diagram,
                               aleph::persistenceIndicatorFunction( pd ),
                               aleph::Envelope()( pd ) } );
        }

        dataSets.push_back( dataSet );
      }
    }
  }

  // Setup distance functor --------------------------------------------

  std::function< double( const PersistenceDiagram&, const PersistenceDiagram&, double ) > functor
    = !useIndicatorFunctionDistance && !useEnvelopeFunctionDistance ?
        useWassersteinDistance ?
          [] ( const PersistenceDiagram& D1, const PersistenceDiagram& D2, double p )
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

  {
    auto name = useEnvelopeFunctionDistance ? "envelope function"
                                            : useIndicatorFunctionDistance
                                              ? "persistence indicator function"
                                              : useWassersteinDistance
                                                ? "Wasserstein"
                                                : "Hausdorff";

    auto type = calculateKernel ? "kernel values" : "distances";

    std::cerr << "* Calculating pairwise " << type << " with " << name << " distance\n";
    std::cerr << "* Calculating pairwise " << type << " with p=" << power << "...";
  }

  std::vector< std::vector<double> > distances;
  distances.resize( dataSets.size(), std::vector<double>( dataSets.size() ) );

  {
    std::size_t n = dataSets.size();
    std::size_t m = dataSets.size() * ( dataSets.size() - 1 ) / 2;

    #pragma omp parallel for
    for( std::size_t k = 0; k < m; k++ )
    {
      auto row = std::size_t( double( n - 2 ) - std::floor( std::sqrt( -8*k + 4*n*(n-1) - 7 ) / 2.0 - 0.5 ) );
      auto col = std::size_t( k + row + 1 - n*(n-1)/2 + (n-row)*( (n-row)-1 ) / 2 );

      double d = 0.0;

      if( useIndicatorFunctionDistance )
        d = distancePIF( dataSets.at(row), dataSets.at(col), minDimension, maxDimension, power, normalize );
      else if( useEnvelopeFunctionDistance )
        d = distanceEnvelopeFunctions( dataSets.at(row), dataSets.at(col), minDimension, maxDimension, power );
      else
        d = persistenceDiagramDistance( dataSets.at(row), dataSets.at(col), minDimension, maxDimension, power, functor );

      if( calculateKernel )
      {
        d = -d;
        if( useExponentialFunction )
          d = std::exp( sigma * d );
      }

      distances[row][col] = d;
      distances[col][row] = d;
    }
  }

  std::cerr << "finished\n";

  std::cerr << "Storing matrix...";

  storeMatrix( distances, std::cout );

  std::cerr << "finished\n";

  std::cerr << "Data sets were processed in the following order:\n";
  for( auto&& dataSet : dataSets )
  {
    if( !dataSet.empty() )
      std::cerr << "  - " << dataSet.front().name << "\n";
  }
}
