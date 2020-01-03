/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  Its purpose is to analyse the persistent homology of connectivity
  matrices, specifically those arising from fMRI data sets. To this
  end, *two* graph filtrations are calculated: one for the positive
  correlations, the other for the negative ones. The resulting data
  will be merged into a single persistence diagram.
*/

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

// Let's make those aliases global because they will be used throughout
// this program.
using DataType           = double;
using VertexType         = unsigned short;
using Simplex            = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex  = aleph::topology::SimplicialComplex<Simplex>;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;
using Point              = typename PersistenceDiagram::Point;

/*
  Class for collecting persistence diagrams for a set of input
  filenames. The class is capable of merging the diagrams that
  correspond to the same filename automatically. Diagrams will
  be merged by taking a union of their points.
*/

class DiagramCollection
{
public:
  template <class InputIterator> DiagramCollection( unsigned numDiagrams, InputIterator begin, InputIterator end )
  {
    // Fill collection with empty diagrams
    for( auto it = begin; it != end; ++it )
      _diagrams[ *it ] = Diagrams( numDiagrams );
  }

private:

  // Slightly verbose, but makes for readable code later on
  using Diagrams = std::vector<PersistenceDiagram>;
  using Value    = Diagrams;
  using Key      = std::string;
  using Map      = std::map<Key, Value>;

  Map _diagrams;
};

void usage()
{
  std::cerr << "Usage: ephemeral [--dimension DIMENSION] [--infinity INF] FILENAMES\n"
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

void processFilename( const std::string& filename,
                      double infinity,
                      bool keepUnpaired,
                      bool verbose,
                      aleph::topology::io::AdjacencyMatrixReader reader )
{
  if( verbose )
    std::cerr << "* Processing " << filename << "...";

  SimplicialComplex K;
  reader( filename, K );

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

    aleph::io::writeJSON( std::cout, diagram, basename, kvs );
  }
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

  std::vector<std::string> filenames;

  for( int i = optind; i < argc; i++ )
    filenames.push_back( argv[i] );

  // The maximum number of diagrams per simplicial complex depends on
  // the maximum expansion dimension and whether we want to keep some
  // unpaired features. This is required for bookkeeping.
  unsigned numDiagrams = keepUnpaired + dimension + 1;

  DiagramCollection diagramCollection(
      numDiagrams,
      filenames.begin(),
      filenames.end()
  );

  aleph::topology::io::AdjacencyMatrixReader reader;
  reader.setIgnoreNaNs();
  reader.setIgnoreZeroWeights();

  // Whew, is there *really* no better way of specifying this strategy
  // here as a qualified name?
  reader.setVertexWeightAssignmentStrategy(
      aleph::topology::io::AdjacencyMatrixReader::VertexWeightAssignmentStrategy::AssignGlobalMinimum
  );

  std::cout << "{\n"
            << "\"diagrams\": [\n";

  for( auto&& filename : filenames )
  {
    processFilename( filename,
                     infinity,
                     keepUnpaired,
                     verbose,
                     reader
    );
  }

  std::cout << "\n"
            << "]\n"
            << "}\n";
}
