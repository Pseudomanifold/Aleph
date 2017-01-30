#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "persistenceDiagrams/IO.hh"
#include "persistenceDiagrams/PersistenceDiagram.hh"
#include "persistenceDiagrams/PersistenceIndicatorFunction.hh"

using DataType           = double;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;
using StepFunction       = aleph::math::StepFunction<double>;

int main( int argc, char** argv )
{
  if( argc <= 3 )
  {
    // TODO: Show usage
    return -1;
  }

  unsigned n = static_cast<unsigned>( std::stoul( argv[1] ) );

  std::vector<std::string> filenames;
  std::vector<PersistenceDiagram> persistenceDiagrams;
  std::vector<StepFunction> persistenceIndicatorFunctions;

  std::set<DataType> domain;

  for( int i = 2; i < argc; i++ )
  {
    filenames.push_back( argv[i] );

    std::cerr << "* Processing '" << argv[i] << "'...";

    PersistenceDiagram persistenceDiagram
        = aleph::load<DataType>( filenames.back() );

    persistenceDiagram.removeDiagonal();
    persistenceDiagram.removeUnpaired();

    persistenceDiagrams.push_back( persistenceDiagram );
    persistenceIndicatorFunctions.push_back( aleph::persistenceIndicatorFunction( persistenceDiagram ) );

    persistenceIndicatorFunctions.back().domain(
      std::inserter( domain, domain.begin() ) );

    std::cerr << "finished\n";
  }

  if( domain.empty() )
    return -1;

  auto min = *domain.begin();
  auto max = *domain.rbegin();

  std::cerr << "* Domain: [" << min << "," << max << "]\n";

  // Prepare bins ------------------------------------------------------

  std::vector<DataType> linbins;
  std::vector<DataType> logbins;

  linbins.reserve( n );
  logbins.reserve( n );

  for( unsigned i = 0; i < n; i++ )
  {
    auto offset = ( max - min ) / (n-1);
    auto value  = min + i * offset;

    linbins.push_back( value );
  }

  std::cerr << "* Linear-spaced bins: ";

  for( auto&& linbin : linbins )
    std::cerr << linbin << " ";

  std::cerr << "\n";

  for( unsigned i = 0; i < n; i++ )
  {
    auto offset = ( std::log10( max ) - std::log10( min ) ) / (n-1);
    auto value  = std::log10( min ) + i * offset;

    logbins.push_back( value );
  }

  std::transform( logbins.begin(), logbins.end(),
                  logbins.begin(),
                  [] ( const DataType x )
                  {
                    return std::pow( DataType(10), x );
                  } );

  std::cerr << "* Log-spaced bins: ";

  for( auto&& logbin : logbins )
    std::cerr << logbin << " ";

  std::cerr << "\n";
}
