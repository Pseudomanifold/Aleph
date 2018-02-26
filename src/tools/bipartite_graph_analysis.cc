#include <aleph/persistenceDiagrams/Norms.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>
#include <aleph/topology/io/BipartiteAdjacencyMatrix.hh>

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

#include <getopt.h>

int main( int argc, char** argv )
{
  bool minimum   = false;
  bool normalize = false;

  {
    static option commandLineOptions[] =
    {
      { "minimum"  , no_argument, nullptr, 'm' },
      { "normalize", no_argument, nullptr, 'n' },
      { nullptr    , 0          , nullptr,  0  }
    };

    int option = 0;
    while( ( option = getopt_long( argc, argv, "mn", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'm':
        minimum = true;
        break;
      case 'n':
        normalize = true;
        break;
      default:
        break;
      }
    }
  }

  using DataType          = double;
  using VertexType        = unsigned short;
  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  DataType minData = std::numeric_limits<DataType>::max();
  DataType maxData = std::numeric_limits<DataType>::lowest();

  // 1. Read simplicial complexes --------------------------------------

  std::vector<SimplicialComplex> simplicialComplexes;
  simplicialComplexes.reserve( static_cast<unsigned>( argc - optind - 1 ) );

  {
    aleph::topology::io::BipartiteAdjacencyMatrixReader reader;

    if( minimum )
      reader.setAssignMinimumVertexWeight();

    for( int i = optind; i < argc; i++ )
    {
      auto filename = std::string( argv[i] );

      std::cerr << "* Processing " << filename << "...";

      SimplicialComplex K;
      reader( filename, K );

      std::cerr << "finished\n";

      // *Always* determine minimum and maximum weights so that we may
      // report them later on. They are only used for normalization in
      // the persistence diagram calculation step.
      for( auto&& s : K )
      {
        minData = std::min( minData, s.data() );
        maxData = std::max( maxData, s.data() );
      }

      simplicialComplexes.emplace_back( K );
    }
  }

  // 2. Calculate persistent homology ----------------------------------

  for( std::size_t i = 0; i < simplicialComplexes.size(); i++ )
  {
    auto&& K      = simplicialComplexes[i];
    auto diagrams = aleph::calculatePersistenceDiagrams( K );
    auto&& D      = diagrams.front();

    D.removeUnpaired();

    if( normalize )
    {
      using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;
      using Point              = typename PersistenceDiagram::Point;

      // Ensures that points are in [-1:+1]
      std::transform( D.begin(), D.end(), D.begin(),
        [&minData, &maxData] ( const Point& p )
        {
          auto x = p.x();
          auto y = p.y();
          x      = 2 * (x - minData) / (maxData - minData) - 1;
          y      = 2 * (y - minData) / (maxData - minData) - 1;

          return Point( x,y );
        }
      );
    }

    std::cout << i << "\t" << aleph::pNorm( D ) << "\n";
  }
}
