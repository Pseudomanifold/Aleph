#include <aleph/math/Bootstrap.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/PersistenceIndicatorFunction.hh>

#include <aleph/persistenceDiagrams/io/Raw.hh>

#include <aleph/utilities/Values.hh>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <cmath>

#include <getopt.h>

using DataType                     = double;
using PersistenceDiagram           = aleph::PersistenceDiagram<DataType>;
using PersistenceIndicatorFunction = decltype( aleph::persistenceIndicatorFunction( PersistenceDiagram() ) );
using Image                        = PersistenceIndicatorFunction::Image;

auto meanCalculation = [] ( auto begin, auto end )
{
  using T  = typename std::iterator_traits<decltype(begin)>::value_type;
  auto sum = std::accumulate( begin, end, T() );

  if( begin == end )
    return std::numeric_limits<T>::quiet_NaN();

  return sum / static_cast<double>( std::distance(begin, end) );
};

unsigned index( unsigned int samples, double alpha )
{
  // This accounts for rounding and works regardless of whether
  // the product samples * alpha is an integer or not. Note the
  // offset of -1. It is required because, say, the 100th value
  // is at index 99 of the vector.
  return static_cast<unsigned>( std::ceil( samples * alpha ) ) - 1;
}

int main( int argc, char** argv )
{
  auto alpha                   = 0.05;
  unsigned numBootstrapSamples = 50;
  bool readStepFunctions       = false;

  {
    static option commandLineOptions[] =
    {
      { "alpha"              , required_argument, nullptr, 'a' },
      { "bootstrap"          , required_argument, nullptr, 'b' },
      { "read-step-functions", no_argument      , nullptr, 's' },
      { nullptr              , 0                , nullptr,  0  }
    };

    int c = 0;
    while( ( c = getopt_long( argc, argv, "a:b:s", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( c )
      {
      case 'a':
        alpha = std::stod( optarg );
        break;

      case 'b':
        numBootstrapSamples = static_cast<unsigned>( std::stoul( optarg ) );
        break;

      case 's':
        readStepFunctions = true;
        break;

      default:
        throw std::runtime_error( "Unknown command-line parameter" );
      }
    }

    aleph::utilities::ensureRange( alpha, 0.0, 1.0 );
    aleph::utilities::ensureLarger( numBootstrapSamples, unsigned(0) );
  }

  // No input files are present, so let's do nothing at all
  if( argc - optind <= 0 )
    return 0;

  std::vector<PersistenceIndicatorFunction> persistenceIndicatorFunctions;
  persistenceIndicatorFunctions.reserve( static_cast<std::size_t>( argc - optind ) );

  for( int i = optind; i < argc; i++ )
  {
    std::cerr << "* Processing '" << argv[i] << "'...";

    if( readStepFunctions )
    {
      PersistenceIndicatorFunction PIF;

      std::ifstream in( argv[i] );
      if( !in )
        throw std::runtime_error( "Unable to load input file" );

      in >> PIF;

      persistenceIndicatorFunctions.emplace_back( PIF );
    }
    else
    {
      auto D = aleph::io::load<DataType>( argv[i] );

      D.removeDiagonal();
      D.removeUnpaired();

      persistenceIndicatorFunctions.emplace_back( aleph::persistenceIndicatorFunction( D ) );
    }

    std::cerr << "finished\n";
  }

  std::vector<PersistenceIndicatorFunction> meanReplicates;
  meanReplicates.reserve( numBootstrapSamples );

  aleph::math::Bootstrap bootstrap;

  bootstrap.makeReplicates( numBootstrapSamples,
                            persistenceIndicatorFunctions.begin(), persistenceIndicatorFunctions.end(),
                            meanCalculation,
                            std::back_inserter( meanReplicates ) );

  auto empiricalMean = meanCalculation( persistenceIndicatorFunctions.begin(), persistenceIndicatorFunctions.end() );

  std::vector<Image> theta;
  theta.reserve( numBootstrapSamples );

  for( auto&& meanReplicate: meanReplicates )
  {
    auto n = persistenceIndicatorFunctions.size();
    auto f = std::sqrt( n ) * ( meanReplicate - empiricalMean );
    f      = f.abs();

    theta.emplace_back( f.sup() );
  }

  std::sort( theta.begin(), theta.end() );

  if( theta.empty() )
    return 0;

  auto quantile = theta.at( index( numBootstrapSamples, 1.0 - alpha ) );
  auto fLower   = empiricalMean - quantile / std::sqrt( persistenceIndicatorFunctions.size() );
  auto fUpper   = empiricalMean + quantile / std::sqrt( persistenceIndicatorFunctions.size() );

  std::ofstream out( "/tmp/Mean_plus_confidence.txt" );

  out << empiricalMean << "\n\n"
      << fUpper        << "\n\n"
      << fLower        << "\n";

}
