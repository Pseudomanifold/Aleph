#include <aleph/math/Bootstrap.hh>
#include <aleph/math/PiecewiseLinearFunction.hh>

#include <aleph/utilities/String.hh>

#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

#include <getopt.h>

#include <cmath>

using DataType = double;
using Function = aleph::math::PiecewiseLinearFunction<DataType>;

Function load( const std::string& filename )
{
  std::ifstream in( filename );
  if( !in )
    throw std::runtime_error( "Unable to open file for reading" );

  std::string line;
  std::vector< std::pair<DataType, DataType> > data;

  while( std::getline( in, line ) )
  {
    line = aleph::utilities::trim( line );

    if( line.empty() || line.front() == '#' )
      continue;

    DataType x = DataType();
    DataType y = DataType();

    std::istringstream converter( line );

    converter >> x >> y;

    if( !converter )
      throw std::runtime_error( "Conversion failed" );

    data.push_back( std::make_pair(x,y) );
  }

  return Function( data.begin(), data.end() );
}

auto meanCalculation = [] ( auto begin, auto end )
{
  using T  = typename std::iterator_traits<decltype(begin)>::value_type;
  auto sum = std::accumulate( begin, end, T() );

  return sum / static_cast<double>( std::distance(begin, end) );
};

int main( int argc, char** argv )
{
  DataType alpha               = 0.05;
  unsigned numBootstrapSamples = 10;

  {
    static option commandLineOptions[] =
    {
      { "alpha", required_argument, nullptr, 'a' },
      { "n"    , required_argument, nullptr, 'n' },
      { nullptr, 0                , nullptr,  0  }
    };

    int option = 0;
    while( ( option = getopt_long( argc, argv, "a:", commandLineOptions, nullptr ) ) != - 1 )
    {
      switch( option )
      {
      case 'a':
        alpha = DataType( std::stod( optarg ) );
        break;
      case 'n':
        numBootstrapSamples = unsigned( std::stoull( optarg ) );
        break;
      default:
        break;
      }
    }
  }

  if( argc - optind <= 1 )
    return -1;

  // Load functions ----------------------------------------------------

  std::cerr << "* Loading functions...";

  std::vector<Function> functions;

  for( int i = optind; i < argc; i++ )
    functions.push_back( load( argv[i] ) );

  std::cerr << "finished\n";

  // Calculate mean ----------------------------------------------------

  std::cerr << "* Calculating empirical mean...";

  // This is the empirical mean that we obtain directly from the input
  // data. We do *not* make any assumptions about its distribution.
  auto empiricalMean = meanCalculation( functions.begin(), functions.end() );

  std::cerr << "finished\n";

  // These are the bootstrap replicates of the mean function. There's
  // one for every bootstrap sample.
  std::vector<Function> meanReplicates;
  meanReplicates.reserve( numBootstrapSamples );

  aleph::math::Bootstrap bootstrap;

  std::cerr << "* Calculating bootstrap replicates (n=" << numBootstrapSamples << ", m=" << functions.size() << ")...";

  bootstrap.makeReplicates( numBootstrapSamples,
                            functions.begin(), functions.end(),
                            meanCalculation,
                            std::back_inserter( meanReplicates ) );

  std::cerr << "finished\n";

  // This contains the population parameter of the corresponding
  // empirical process, viz. the *supremum* of the difference in
  // empirical mean and bootstrapped mean.
  std::vector<DataType> theta;
  theta.reserve( numBootstrapSamples );

  std::cerr << "* Calculating confidence band information...";

  for( auto&& meanReplicate : meanReplicates )
  {
    auto n = functions.size();
    auto f = std::sqrt( n ) * ( meanReplicate - empiricalMean );
    f      = f.abs();

    theta.push_back( f.sup() );
  }

  std::sort( theta.begin(), theta.end() );

  std::cerr << "finished\n";

  auto quantile = theta.at( bootstrap.index( numBootstrapSamples, alpha / 2 ) );
  auto fLower   = empiricalMean - quantile / std::sqrt( functions.size() );
  auto fUpper   = empiricalMean + quantile / std::sqrt( functions.size() );

  std::cout << empiricalMean << "\n\n"
            << fUpper        << "\n\n"
            << fLower        << "\n";
}
