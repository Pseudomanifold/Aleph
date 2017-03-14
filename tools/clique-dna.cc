#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <getopt.h>

#include "persistenceDiagrams/IO.hh"
#include "persistenceDiagrams/PersistenceDiagram.hh"
#include "persistenceDiagrams/PersistenceIndicatorFunction.hh"

using DataType           = double;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;
using StepFunction       = aleph::math::StepFunction<double>;

/*
  Prints a vector in a simple matrix-like format. Given a row index, all
  of the entries are considered to be the columns of the matrix:

    Input:   {4,5,6}, row = 23

    Output:  23 0 4
             23 1 5
             23 2 6
             <empty line>

  This output format is flexible and permits direct usage in other tools
  such as TikZ or pgfplots.
*/

template <class Map> void print( std::ostream& o, const Map& m, unsigned row )
{
  unsigned column = 0;
  for( auto it = m.begin(); it != m.end(); ++it )
    o << column++ << "\t" << row << "\t" << *it << "\n";

  o << "\n";
}

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "min-k", required_argument, nullptr, 'k' },
    { "max-k", required_argument, nullptr, 'K' },
    { nullptr, 0                , nullptr,  0  }
  };

  int option = 0;

  unsigned minK = 0;
  unsigned maxK = 0;

  while( ( option = getopt_long( argc, argv, "k:K:", commandLineOptions, nullptr ) ) != -1 )
  {
    switch( option )
    {
    case 'k':
      minK = unsigned( std::stoul( optarg ) );
      break;
    case 'K':
      maxK = unsigned( std::stoul( optarg ) );
      break;
    default:
      break;
    }
  }

  if( argc - optind < 2 )
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

  unsigned index = 0;
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

      if( logbin < loghist.size() )
        loghist.at(logbin) += value;
    }

    print( linout, linhist, index );
    print( logout, loghist, index );

    ++index;
  }
}
