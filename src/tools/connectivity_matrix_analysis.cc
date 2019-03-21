/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  Its purpose is to analyse the persistent homology of connectivity
  matrices. The tool bears some semblance to the *network analysis*
  tools, but focuses specifically on data sets whose weights are an
  interpretable correlation measure.
*/

#include <aleph/persistenceDiagrams/Entropy.hh>
#include <aleph/persistenceDiagrams/Norms.hh>

#include <aleph/persistenceDiagrams/io/JSON.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/io/AdjacencyMatrix.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/utilities/Filesystem.hh>

#include <cmath>

#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <getopt.h>

void usage()
{
  std::cerr << "Usage: connectivity_matrix_analysis [--dimension DIMENSION] [--infinity INF] FILENAMES\n"
            << "\n"
            << "Analyses a set of connectivity matrices. The matrices are optionally\n"
            << "expanded to a pre-defined dimension. By default, only information of\n"
            << "the zeroth persistent homology group will be shown.\n"
            << "\n"
            << "The value INF will be used to replace infinite values in the diagram\n"
            << "in order to facilitate the subsequent analysis.\n"
            << "\n"
            << "Flags:\n"
            << "  -k: keep & report unpaired simplices (infinite values)\n"
            << "  -v: verbose output\n"
            << "\n";
}

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "dimension"     , required_argument, nullptr, 'd' },
    { "infinity"      , required_argument, nullptr, 'i' },
    { "keep-unpaired" , no_argument      , nullptr, 'k' },
    { "verbose"       , no_argument      , nullptr, 'v' },
    { nullptr         , 0                , nullptr,  0  }
  };

  unsigned dimension = 0;
  double infinity    = std::numeric_limits<double>::infinity();
  bool keepUnpaired  = false;
  bool verbose       = false;

  {
    int option = 0;

    while( ( option = getopt_long( argc, argv, "d:i:kv", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'd':
        dimension = static_cast<unsigned>( std::stoul(optarg) );
        break;
      case 'i':
        infinity = static_cast<double>( std::stod(optarg) );
        break;
      case 'k':
        keepUnpaired = true;
        break;
      }
    }
  }

  if( (argc - optind) < 1 )
  {
    usage();
    return -1;
  }

  using DataType           = double;
  using VertexType         = unsigned short;
  using Simplex            = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex  = aleph::topology::SimplicialComplex<Simplex>;
  using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;
  using Point              = typename PersistenceDiagram::Point;

  std::vector<std::string> filenames;

  for( int i = optind; i < argc; i++ )
    filenames.push_back( argv[i] );

  aleph::topology::io::AdjacencyMatrixReader reader;
  reader.setIgnoreNaNs();
  reader.setIgnoreZeroWeights();

  // Whew, is there *really* no better way of specifying this strategy
  // here as a qualified name?
  reader.setVertexWeightAssignmentStrategy(
      aleph::topology::io::AdjacencyMatrixReader::VertexWeightAssignmentStrategy::AssignZero
  );

  std::cout << "{\n"
            << "\"diagrams\": [\n";

  for( auto&& filename : filenames )
  {
    if( verbose )
      std::cerr << "* Processing " << filename << "...";

    SimplicialComplex K;
    reader( filename, K,
      [] ( DataType maxWeight, DataType /* minWeight */, DataType weight )
      {
        // Transform the weight into a *distance* by negating it; this
        // ignores all other scaling mechanisms applied to the data.
        return maxWeight - weight;
      }
    );

    K.sort();

    bool dualize                    = true;
    bool includeAllUnpairedCreators = keepUnpaired;

    auto diagrams
      = aleph::calculatePersistenceDiagrams( K,
                                             dualize,
                                             includeAllUnpairedCreators );

    if( verbose )
      std::cerr << "finished\n";

    auto basename
      = aleph::utilities::basename( filename );

    for( auto&& diagram : diagrams )
    {
      if( std::isfinite( infinity ) )
      {
        std::transform( diagram.begin(), diagram.end(), diagram.begin(),
            [&infinity] ( const Point& p )
            {
              if( p.isUnpaired() )
                return Point( p.x(), infinity );
              else
                return Point( p.x(), p.y() );
            }
        );
      }

      // Stores additional data about each persistence diagram in order
      // to make it easier to keep track of information.
      std::map<std::string, std::string> kvs;

      kvs["total_persistence_1"] = std::to_string( aleph::totalPersistence( diagram, 1.0 ) );
      kvs["total_persistence_2"] = std::to_string( aleph::totalPersistence( diagram, 2.0 ) );

      kvs["persistent_entropy"]  = std::to_string( aleph::persistentEntropy( diagram ) );

      aleph::io::writeJSON( std::cout, diagram, basename, kvs );
    }
  }

  std::cout << "\n"
            << "]\n"
            << "}\n";
}
