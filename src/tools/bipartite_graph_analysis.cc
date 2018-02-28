#include <aleph/math/SymmetricMatrix.hh>

#include <aleph/persistenceDiagrams/Distances.hh>
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
  bool absolute              = false;
  bool minimum               = false;
  bool normalize             = false;
  bool calculateDiagrams     = false;
  bool calculateTrajectories = false;

  {
    static option commandLineOptions[] =
    {
      { "absolute"            , no_argument, nullptr, 'a' },
      { "minimum"             , no_argument, nullptr, 'm' },
      { "normalize"           , no_argument, nullptr, 'n' },
      { "persistence-diagrams", no_argument, nullptr, 'p' },
      { "trajectories"        , no_argument, nullptr, 't' },
      { nullptr               , 0          , nullptr,  0  }
    };

    int option = 0;
    while( ( option = getopt_long( argc, argv, "amnpt", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'a':
        absolute = true;
        break;
      case 'm':
        minimum = true;
        break;
      case 'n':
        normalize = true;
        break;
      case 'p':
        calculateDiagrams = true;
        break;
      case 't':
        calculateTrajectories = true;
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

    if( absolute )
      reader.setUseAbsoluteValues();

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

  using Matrix = aleph::math::SymmetricMatrix<double>;

  // Stores the distance matrix for the trajectories of persistence
  // diagrams. This will only be used if the client set the correct
  // flag.
  Matrix trajectoryDistances;
  if( calculateTrajectories )
    trajectoryDistances = Matrix( simplicialComplexes.size() );

  using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;
  using Point              = typename PersistenceDiagram::Point;

  // Stores the zeroth persistence diagram for calculating trajectories
  // later on. This may need to be extended in order to handle diagrams
  // with higher-dimensional features.
  std::vector<PersistenceDiagram> trajectoryDiagrams;
  if( calculateTrajectories )
    trajectoryDiagrams.reserve( simplicialComplexes.size() );

  for( std::size_t i = 0; i < simplicialComplexes.size(); i++ )
  {
    auto&& K      = simplicialComplexes[i];
    auto diagrams = aleph::calculatePersistenceDiagrams( K );
    auto&& D      = diagrams.front();

    D.removeDiagonal();
    D.removeUnpaired();

    if( normalize )
    {
      // Ensures that points are in [-1:+1] or [0:1] if absolute values
      // have been selected by the client.
      std::transform( D.begin(), D.end(), D.begin(),
        [&absolute, &minData, &maxData] ( const Point& p )
        {
          auto x = p.x();
          auto y = p.y();

          if( absolute )
          {
            x = (x - minData) / (maxData - minData);
            y = (y - minData) / (maxData - minData);
          }
          else
          {
            x = 2 * (x - minData) / (maxData - minData) - 1;
            y = 2 * (y - minData) / (maxData - minData) - 1;
          }

          return Point( x,y );
        }
      );
    }

    // Determine mode of operation -------------------------------------
    //
    // Several modes of operation exist for this program. They can be
    // set using the flags specified above. At present, the following
    // operations are possible:
    //
    // - Calculate persistence diagrams
    // - Calculate persistence diagram trajectories
    // - Calculate 2-norm of the persistence diagrams

    if( calculateDiagrams )
      std::cout << D << "\n\n";
    else if( calculateTrajectories )
      trajectoryDiagrams.push_back( D );
    else
      std::cout << i << "\t" << aleph::pNorm( D ) << "\n";
  }

  // Need to calculate the trajectories afterwards because they require
  // building a database of persistence diagrams.
  if( calculateTrajectories )
  {
    for( std::size_t i = 0; i < trajectoryDiagrams.size(); i++ )
    {
      auto&& Di = trajectoryDiagrams[i];

      for( std::size_t j = i+1; j < trajectoryDiagrams.size(); j++ )
      {
        auto&& Dj = trajectoryDiagrams[j];
        auto dist = aleph::distances::hausdorffDistance(
          Di, Dj
        );

        trajectoryDistances(i,j) = dist;
      }
    }

    // FIXME: replace with proper layout
    std::cout << trajectoryDistances;
  }
}
