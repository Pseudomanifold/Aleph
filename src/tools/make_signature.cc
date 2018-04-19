#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>
#include <aleph/persistenceDiagrams/io/Raw.hh>

#include <algorithm>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <getopt.h>

using DataType           = double;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;
using Point              = typename PersistenceDiagram::Point;

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

int main( int argc, char** argv)
{
  bool normalize = false;

  {
    static option commandLineOptions[] =
    {
      { "normalize", no_argument, nullptr, 'n' },
      { nullptr    , 0          , nullptr,  0  }
    };

    int option = 0;

    while( ( option = getopt_long( argc, argv, "n", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
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

    if( normalize )
      normalizeDiagram( diagram );

    diagrams.emplace_back( diagram );

    std::cerr << "finished\n";
  }
}
