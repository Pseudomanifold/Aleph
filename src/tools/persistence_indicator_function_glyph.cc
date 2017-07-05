/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  Given a set of persistence diagrams, the tool calculates a histogram
  glyph. The glyph uses persistence indicator functions, a summarizing
  function of a persistence diagram.

  This tool follows the publication:

    Clique Community Persistence: A Topological Visual Analysis Approach for Complex Networks
    Bastian Rieck, Ulderico Fugacci, Jonas Lukasczyk, Heike Leitte
    Submitted to IEEE Vis 2017

  TODO: Link to documentation
*/

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <cmath>

#include <getopt.h>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/PersistenceIndicatorFunction.hh>

#include <aleph/persistenceDiagrams/io/Raw.hh>

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
    { "max-k", required_argument, nullptr, 'K' },
    { nullptr, 0                , nullptr,  0  }
  };

  int option    = 0;
  unsigned maxK = 0;

  while( ( option = getopt_long( argc, argv, "K:", commandLineOptions, nullptr ) ) != -1 )
  {
    switch( option )
    {
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

  unsigned n = static_cast<unsigned>( std::stoul( argv[optind++] ) );

  std::vector<std::string> filenames;
  std::vector<PersistenceDiagram> persistenceDiagrams;
  std::vector<StepFunction> persistenceIndicatorFunctions;

  std::set<DataType> domain;

  for( int i = optind; i < argc; i++ )
  {
    filenames.push_back( argv[i] );

    std::cerr << "* Processing '" << argv[i] << "'...";

    PersistenceDiagram persistenceDiagram
        = aleph::io::load<DataType>( filenames.back() );

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
  auto max = std::nextafter( *domain.rbegin(), std::numeric_limits<DataType>::max() );

  std::cerr << "* Domain: [" << min << "," << max << "]\n";

  // Prepare bins ------------------------------------------------------

  std::vector<DataType> linbins;
  std::vector<DataType> logbins;

  linbins.reserve( n );
  logbins.reserve( n );

  for( unsigned i = 0; i < n; i++ )
  {
    auto offset = ( max - min ) / n;
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

  for( auto it = linbins.begin(); it != linbins.end(); ++it )
  {
    if( it != linbins.begin() )
    {
      auto prev = std::prev( it );
      domain.insert( ( *prev + *it ) / 2.0 );
    }
  }

  for( auto it = logbins.begin(); it != logbins.end(); ++it )
  {
    if( it != logbins.begin() && std::isfinite( *it ) )
    {
      auto prev = std::prev( it );
      domain.insert( ( *prev + *it ) / 2.0 );
    }
  }

  // Replace domain ---------------------------------------------------

  {
    std::set<DataType> newDomain;

    for( auto&& x : domain )
    {
      if( x != *domain.rbegin() )
        newDomain.insert( std::nextafter( x,  std::numeric_limits<DataType>::max() ) );
      if( x != *domain.begin() )
        newDomain.insert( std::nextafter( x, -std::numeric_limits<DataType>::max() ) );
    }

    domain.swap( newDomain );
  }

  // Prepare histogram calculation -------------------------------------

  auto valueToLinIndex = [&min, &max, &linbins, &n] ( DataType value )
  {
    auto offset = ( max - min ) / n;
    return static_cast<std::size_t>( ( value - min ) / offset );
  };

  auto valueToLogIndex = [&min, &max, &logbins, &n] ( DataType value )
  {
    auto offset = ( std::log10( max ) - std::log10( min ) ) / n;
    return static_cast<std::size_t>( ( std::log10( value ) - std::log10( min ) ) / offset );
  };

  std::ofstream linout( "/tmp/Persistence_indicator_function_glyph_" + std::to_string( n ) + "_lin.txt" );
  std::ofstream logout( "/tmp/Persistence_indicator_function_glyph_" + std::to_string( n ) + "_log.txt" );

  unsigned index = 0;
  for( auto&& pif : persistenceIndicatorFunctions )
  {
    std::vector<DataType> linhist(n);
    std::vector<DataType> loghist(n);

    for( auto&& x : domain )
    {
      auto value  = pif(x);
      auto linbin = valueToLinIndex(x);
      auto logbin = valueToLogIndex(x);

      linhist.at(linbin) = std::max( linhist.at(linbin), value );

      if( logbin < loghist.size() )
        loghist.at(logbin) = std::max( loghist.at(logbin), value );
    }

    print( linout, linhist, index );
    print( logout, loghist, index );

    ++index;
  }

  // Extend the output with sufficiently many empty histograms. This
  // ensures that the output has the same dimensions.
  for( unsigned i = unsigned( persistenceIndicatorFunctions.size() ); i < maxK; i++ )
  {
    std::vector<DataType> hist(n);

    print( linout, hist, i );
    print( logout, hist, i );
  }
}
