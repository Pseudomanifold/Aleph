#include <aleph/math/Bootstrap.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/PersistenceIndicatorFunction.hh>

#include <aleph/persistenceDiagrams/io/Raw.hh>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <cmath>

using DataType                     = double;
using PersistenceDiagram           = aleph::PersistenceDiagram<DataType>;
using PersistenceIndicatorFunction = decltype( aleph::persistenceIndicatorFunction( PersistenceDiagram() ) );
using Image                        = PersistenceIndicatorFunction::Image;

auto meanCalculation = [] ( auto begin, auto end )
{
  using T  = typename std::iterator_traits<decltype(begin)>::value_type;
  auto sum = std::accumulate( begin, end, T() );

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
  std::vector<PersistenceIndicatorFunction> persistenceIndicatorFunctions;
  persistenceIndicatorFunctions.reserve( static_cast<std::size_t>( argc - 1 ) );

  for( int i = 1; i < argc; i++ )
  {
    std::cerr << "* Processing '" << argv[i] << "'...";

    auto D = aleph::io::load<DataType>( argv[i] );

    D.removeDiagonal();
    D.removeUnpaired();

    persistenceIndicatorFunctions.emplace_back( aleph::persistenceIndicatorFunction( D ) );

    std::cerr << "finished\n";
  }

  // TODO: make configurable
  unsigned numBootstrapSamples = 50;

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

  auto alpha    = 0.05;
  auto quantile = theta.at( index( numBootstrapSamples, alpha / 2 ) );
  auto fLower   = empiricalMean - quantile / std::sqrt( persistenceIndicatorFunctions.size() );
  auto fUpper   = empiricalMean + quantile / std::sqrt( persistenceIndicatorFunctions.size() );

  std::ofstream out( "/tmp/Mean_plus_confidence.txt" );

  out << empiricalMean << "\n\n"
      << fUpper        << "\n\n"
      << fLower        << "\n";

}
