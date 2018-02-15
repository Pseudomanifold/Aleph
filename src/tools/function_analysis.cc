#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/topology/io/Function.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/utilities/String.hh>

#include <fstream>
#include <istream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <cmath>

#include <getopt.h>

using DataType           = double;
using VertexType         = unsigned;
using Simplex            = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex  = aleph::topology::SimplicialComplex<Simplex>;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

/**
  Auxiliary function for extracting the minimum and maximum data value
  of a given simplicial complex.
*/

std::pair<DataType, DataType> minmaxData( const SimplicialComplex& K )
{
  DataType min = std::numeric_limits<DataType>::max();
  DataType max = std::numeric_limits<DataType>::lowest();

  for( auto&& s : K )
  {
    auto d = s.data();
    min    = std::min( min, d );
    max    = std::max( max, d );
  }

  return std::make_pair( min, max );
}

std::vector<SimplicialComplex> readData( std::istream& in, bool useSublevelSetFiltration )
{
  std::vector<SimplicialComplex> complexes;
  std::string line;

  while( std::getline( in, line ) )
  {
    auto tokens
      = aleph::utilities::split(
          line,
          std::string( "[:;,[:space:]]+" )
    );

    std::vector<DataType> values;
    values.reserve( tokens.size() );

    for( auto&& token : tokens )
    {
      bool success = false;
      auto value   = aleph::utilities::convert<DataType>( token, success );

      if( !success )
        throw std::runtime_error( "Unable to convert token to expected data type" );

      values.emplace_back( value );
    }

    complexes.push_back(
      aleph::topology::io::loadFunction<SimplicialComplex>(
        values.begin(), values.end(),
        [&useSublevelSetFiltration] ( DataType x, DataType y )
        {
          if( useSublevelSetFiltration )
            return std::max(x,y);
          else
            return std::min(x,y);
        }
      )
    );

    auto&& K = complexes.back();

    // Establish filtration order of the simplicial complex. For the
    // sublevel set filtration, regular sorting is sufficient, while
    // for the superlevel set filtration, the comparison functor has
    // to be swapped out.
    if( useSublevelSetFiltration )
    {
      K.sort(
        aleph::topology::filtrations::Data< Simplex, std::less<DataType> >()
      );
    }
    else
    {
      K.sort(
        aleph::topology::filtrations::Data< Simplex, std::greater<DataType> >()
      );
    }
  }

  return complexes;
}

void usage()
{
}

int main( int argc, char** argv )
{
  // Options parsing ---------------------------------------------------
  //
  // By default, a sublevel set filtration is being calculated for the
  // input data set. One or more input data sets may be specified at a
  // time. Using '-' indicates that input should be read from `stdin`.

  bool useSublevelSetFiltration = true;

  {
    static option commandLineOptions[] =
    {
      { "sublevels"  , no_argument, nullptr, 's' },
      { "superlevels", no_argument, nullptr, 'S' },
      { nullptr      , 0          , nullptr,  0  }
    };

    int option = 0;
    while( ( option = getopt_long( argc, argv, "sS", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 's':
        useSublevelSetFiltration = true;
        break;

      case 'S':
        useSublevelSetFiltration = false;
        break;

      default:
        break;
      }
    }
  }

  if( ( argc - optind ) < 1 )
  {
    usage();
    return -1;
  }

  std::vector<SimplicialComplex> complexes;

  for( int i = optind; i < argc; i++ )
  {
    std::string filename = argv[i];
    std::cerr << "* Reading '" << filename << "'...";

    std::istream* in = &std::cin;
    std::ifstream fin;

    if( filename != "-" and not filename.empty() )
    {
      fin.open( filename );
      in = &fin;
    }

    auto localComplexes
      = readData( *in, useSublevelSetFiltration );

    std::cerr << "finished\n";

    complexes.insert( complexes.end(),
      localComplexes.begin(), localComplexes.end()
    );
  }

  std::cerr << "* Read " << complexes.size() << " simplicial complexes\n";

  // Persistent homology calculation -----------------------------------
  //
  // Calculate the zero-dimensional persistent homology of every stored
  // complex.

  std::cerr << "* Calculating persistent homology...";

  for( auto&& K : complexes )
  {
    auto diagrams = aleph::calculatePersistenceDiagrams( K );
    auto minmax   = minmaxData( K );

    if( diagrams.size() != 1 )
      throw std::runtime_error( "Unexpected number of persistence diagrams" );

    auto&& D    = diagrams.front();
    using Point = typename PersistenceDiagram::Point;

    if( D.betti() != 1 )
      throw std::runtime_error( "Unexpected Betti number" );

    std::transform( D.begin(), D.end(), D.begin(),
      [&minmax, &useSublevelSetFiltration] ( const Point& p )
      {
        if( !std::isfinite( p.y() ) )
        {
          // Use the *maximum* weight for the sublevel set filtration so
          // that all points are *above* the diagonal, and vice versa in
          // case of the superlevel set filtration.
          auto y = useSublevelSetFiltration
            ? minmax.second
            : minmax.first;

          return Point( p.x(), y );
        }

        // Just copy the original point; this is not highly efficient
        // but the amount of data should not be too large
        else
          return Point( p );
      }
    );

    std::cout << D << std::endl;
  }

  std::cerr << "finished\n";
}
