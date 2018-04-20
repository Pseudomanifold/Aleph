#include <aleph/geometry/distances/Infinity.hh>

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/io/Raw.hh>

#include <cmath>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <getopt.h>

using DataType           = double;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;
using Point              = typename PersistenceDiagram::Point;
using InfinityDistance   = aleph::geometry::distances::InfinityDistance<DataType>;

void filterDiagram( PersistenceDiagram& D, DataType threshold, bool lower = true )
{
  D.erase(
    std::remove_if( D.begin(), D.end(),
      [&threshold, &lower] ( const Point& p )
      {
        // Note that the *absolute* persistence value should be used
        // here because it will always work, regardless of the state
        // of the corresponding filtration that was employed for the
        // persistence diagram calculation.
        if( lower )
          return std::abs( p.persistence() ) < threshold;
        else
          return std::abs( p.persistence() ) > threshold;
      }
    ),
    D.end()
  );
}

void normalizeDiagram( PersistenceDiagram& D )
{
  auto min = std::numeric_limits<DataType>::max();
  auto max = std::numeric_limits<DataType>::lowest();

  if( D.betti() != 0 )
    throw std::runtime_error( "Normalization not yet implemented for unpaired points" );

  for( auto&& p : D )
  {
    min = std::min( min, std::min( p.x(), p.y() ) );
    max = std::max( max, std::max( p.x(), p.y() ) );
  }

  // Silently ignore invalid ranges of persistence diagram points
  // because they do not influence the results.
  if( min == max )
    return;

  std::transform( D.begin(), D.end(), D.begin(),
    [&min, &max] ( const Point& p )
    {
      auto x = (p.x() - min) / (max-min);
      auto y = (p.y() - min) / (max-min);

      return Point(x, y);
    }
  );
}

std::vector<DataType> makeSignature( const PersistenceDiagram& D )
{
  InfinityDistance distance;

  // All pairwise distances with potential repetitions for points that
  // are well outside the influence radius of other points.
  std::vector<DataType> distances;
  distances.reserve( ( D.size() * D.size() - 1 ) / 2 );

  for( auto it1 = D.begin(); it1 != D.end(); ++it1 )
  {
    for( auto it2 = std::next( it1 ); it2 != D.end(); ++it2 )
    {
      auto dxy = distance( *it1, *it2 );
      auto dx  = std::abs( it1->persistence() );
      auto dy  = std::abs( it2->persistence() );

      distances.emplace_back( std::min( dxy, std::min( dx, dy ) ) );
    }
  }

  std::sort( distances.begin(), distances.end(),
             std::greater<DataType>() );

  return distances;
}

template <class InputIterator> void printSignature(
  InputIterator begin,
  InputIterator end,
  std::ostream& out )
{
  for( auto it = begin; it != end; ++it )
  {
    if( it != begin )
      out << " ";

    out << *it;
  }

  out << "\n";
}

int main( int argc, char** argv)
{
  bool normalize     = false;
  unsigned keep      = 0;
  DataType threshold = DataType();

  {
    static option commandLineOptions[] =
    {
      { "filter"   , required_argument, nullptr, 'f' },
      { "keep"     , required_argument, nullptr, 'k' },
      { "normalize", no_argument      , nullptr, 'n' },
      { nullptr    , 0                , nullptr,  0  }
    };

    int option = 0;

    while( ( option = getopt_long( argc, argv, "f:k:n", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'f':
        threshold = static_cast<DataType>( std::stod( optarg ) );
        break;
      case 'k':
        keep = static_cast<unsigned>( std::stoul( optarg ) );
        break;
      case 'n':
        normalize = true;
        break;
      default:
        break;
      }
    }
  }

  std::vector<PersistenceDiagram> diagrams;
  diagrams.reserve( static_cast<std::size_t>( argc - optind ) );

  for( int i = optind; i < argc; i++ )
  {
    std::cerr << "* Loading '" << argv[i] << "'...";

    auto diagram
      = aleph::io::load<DataType>( argv[i] );

    std::cerr << "finished\n"
              << "* Loaded diagram with " << diagram.size() << " points\n";

    if( threshold != DataType() )
    {
      std::cerr << "* Filtering diagram...";

      filterDiagram( diagram, threshold );

      std::cerr << "finished\n"
                << "* Filtered diagram contains " << diagram.size() << " points\n";
    }

    if( normalize )
    {
      std::cerr << "* Normalizing diagram...";

      normalizeDiagram( diagram );

      std::cerr << "finished\n";
    }

    diagrams.emplace_back( diagram );

  }

  for( auto&& diagram : diagrams )
  {
    auto signature = makeSignature( diagram );

    if( keep != 0 )
      signature.resize( keep );

    printSignature( signature.begin(), signature.end(), std::cout );
  }
}
