#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "persistenceDiagrams/IO.hh"
#include "persistenceDiagrams/PersistenceDiagram.hh"
#include "persistenceDiagrams/PersistenceIndicatorFunction.hh"

using DataType           = double;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;
using StepFunction       = aleph::math::StepFunction<double>;

template <class Map> void print( std::ostream& o, const Map& m )
{
  for( auto it = m.begin(); it != m.end(); ++it )
  {
    if( it != m.begin() )
      o << "\t";

    o << *it;
  }

  o << "\n";
}

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

  // Prepare histogram calculation -------------------------------------

  auto valueToLinIndex = [&min, &max, &linbins, &n] ( DataType value )
  {
    auto offset = ( max - min ) / (n-1);
    return static_cast<std::size_t>( ( value - min ) / offset );
  };

  auto valueToLogIndex = [&min, &max, &logbins, &n] ( DataType value )
  {
    auto offset = ( std::log10( max ) - std::log10( min ) ) / (n-1);
    return static_cast<std::size_t>( ( std::log10( value ) - std::log10( min ) ) / offset );
  };

  std::ofstream linout( "/tmp/DNA_" + std::to_string( n ) + "_lin.txt" );
  std::ofstream logout( "/tmp/DNA_" + std::to_string( n ) + "_log.txt" );

  for( auto&& pif : persistenceIndicatorFunctions )
  {
    std::vector<DataType> linhist(n);
    std::vector<DataType> loghist(n);

    std::set<DataType> domain;
    pif.domain( std::inserter( domain, domain.begin() ) );

    for( auto&& x : domain )
    {
      auto value  = pif(x);
      auto linbin = valueToLinIndex(x);
      auto logbin = valueToLogIndex(x);

      linhist.at(linbin) += value;
      loghist.at(logbin) += value;
    }

    print( linout, linhist );
    print( logout, loghist );
  }
}
