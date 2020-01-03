/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  Its purpose is to analyse the persistent homology of connectivity
  matrices, specifically those arising from fMRI data sets. To this
  end, *two* graph filtrations are calculated: one for the positive
  correlations, the other for the negative ones. The resulting data
  will be merged into a single persistence diagram.
*/

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/topology/io/AdjacencyMatrix.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/utilities/Filesystem.hh>

#include <cassert>
#include <cmath>

#include <algorithm>
#include <iostream>
#include <fstream>
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

  // Slightly verbose, but makes for readable code later on
  using Diagrams = std::vector<PersistenceDiagram>;
  using Value    = Diagrams;
  using Key      = std::string;
  using Map      = std::map<Key, Value>;

  template <class InputIterator> DiagramCollection( unsigned numDiagrams, InputIterator begin, InputIterator end )
  {
    // Fill collection with empty diagrams
    for( auto it = begin; it != end; ++it )
      _diagrams[ *it ] = Diagrams( numDiagrams );
  }

  template <class InputIterator> void update( const Key& key,
                                              InputIterator begin,
                                              InputIterator end )
  {
    // Let's not check for existence here but gamble a little bit ;-)
    auto&& diagrams = _diagrams.at( key );
    unsigned index  = 0;

    // Again, no checking bounds because I am a little bit lazy.
    for( auto it = begin; it != end; ++it )
      diagrams.at( index++ ).merge( *it );
  }

  Diagrams& operator[]( const Key& key )
  {
    return _diagrams[ key ];
  }

private:

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

std::vector<PersistenceDiagram> processFilename( const std::string& filename,
                                                 double infinity,
                                                 bool keepUnpaired,
                                                 bool verbose,
                                                 bool reverse,
                                                 bool distance,
                                                 unsigned numDiagrams,
                                                 aleph::topology::io::AdjacencyMatrixReader reader )
{
  if( verbose )
    std::cerr << "* Processing " << filename << "...";

  SimplicialComplex K;

  if( distance )
  {
    reader( filename, K,
      [] ( DataType /* maxWeight */, DataType /* minWeight */, DataType weight )
      {
        // Transform the weight into a *distance* by negating it; this
        // ignores all other scaling mechanisms applied to the data.
        return 1.0 - weight;
      }
    );
  }
  else
    reader( filename, K );

  // Setting both of them would be invalid
  assert( !(distance && reverse) );

  if( reverse )
  {
    K.sort(
       aleph::topology::filtrations::Data< Simplex, std::less<DataType> >()
     );
  }
  else
    K.sort( aleph::topology::filtrations::Data<Simplex>() );

  bool dualize                    = true;
  bool includeAllUnpairedCreators = keepUnpaired;

  auto diagrams
    = aleph::calculatePersistenceDiagrams( K,
                                           dualize,
                                           includeAllUnpairedCreators );

  diagrams.resize( numDiagrams );

  // Ensures that non-empty diagrams follow the indexing of the vector.
  // For example, if we have data for dimensions 0 and 1, the diagrams
  // should be stored at index 0 and 1, respectively. This is a sanity
  // check that will fail if the data behaves weirdly.
  for( std::size_t i = 0; i < diagrams.size(); i++ )
  {
    auto&& D = diagrams[i];

    if( !D.empty() )
      assert( D.dimension() == i );
  }

  if( verbose )
    std::cerr << "finished\n";

  auto basename
    = aleph::utilities::basename( filename );

  // Negate any finite values supplied by the user to ensure symmetry of
  // the reverse filtration.
  if( reverse && std::isfinite( infinity ) )
    infinity = -infinity;

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
  }

  return diagrams;
}

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "dimension"     , required_argument, nullptr, 'd' },
    { "infinity"      , required_argument, nullptr, 'i' },
    { "keep-unpaired" , no_argument      , nullptr, 'k' },
    { "verbose"       , no_argument      , nullptr, 'v' },
    { "distance"      , no_argument      , nullptr, 'D' },
    { nullptr         , 0                , nullptr,  0  }
  };

  unsigned dimension = 0;
  double infinity    = std::numeric_limits<double>::infinity();
  bool keepUnpaired  = false;
  bool verbose       = false;
  bool distance      = false;

  {
    int option = 0;

    while( ( option = getopt_long( argc, argv, "d:i:kvD", commandLineOptions, nullptr ) ) != -1 )
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
      case 'D':
        distance = true;
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

  // No distance calculations are desired; calculate dual filtration and
  // store them.
  if( !distance )
  {
    // Ascending filtration ----------------------------------------------
    //
    // This filtration goes from *negatively* correlated features of the
    // graphs to positively correlated ones.

    bool reverse = false;

    // Whew, is there *really* no better way of specifying this strategy
    // here as a qualified name?
    reader.setVertexWeightAssignmentStrategy(
        aleph::topology::io::AdjacencyMatrixReader::VertexWeightAssignmentStrategy::AssignGlobalMinimum
    );

    for( auto&& filename : filenames )
    {
      auto diagrams = processFilename( filename,
                                       infinity,
                                       keepUnpaired,
                                       verbose,
                                       reverse,
                                       distance,
                                       numDiagrams,
                                       reader
      );

      diagramCollection.update( filename, diagrams.begin(), diagrams.end() );
    }

    // Descending filtration ---------------------------------------------
    //
    // This filtration goes from *positively* correlated features of the
    // graphs to negatively correlated ones.

    reverse = true;

    // Whew, is there *really* no better way of specifying this strategy
    // here as a qualified name?
    reader.setVertexWeightAssignmentStrategy(
        aleph::topology::io::AdjacencyMatrixReader::VertexWeightAssignmentStrategy::AssignGlobalMaximum
    );

    for( auto&& filename : filenames )
    {
      auto diagrams = processFilename( filename,
                                       infinity,
                                       keepUnpaired,
                                       verbose,
                                       reverse,
                                       distance,
                                       numDiagrams,
                                       reader
      );

      diagramCollection.update( filename, diagrams.begin(), diagrams.end() );
    }
  }

  // Distance calculations are desired, rephrase the expansion and
  // creation of simplicial complexes accordingly.
  else
  {
    reader.setVertexWeightAssignmentStrategy(
        aleph::topology::io::AdjacencyMatrixReader::VertexWeightAssignmentStrategy::AssignZero
    );

    for( auto&& filename : filenames )
    {
      auto diagrams = processFilename( filename,
                                       infinity,
                                       keepUnpaired,
                                       verbose,
                                       false, // no reverse filtration
                                       true,  // distance
                                       numDiagrams,
                                       reader
      );

      diagramCollection.update( filename, diagrams.begin(), diagrams.end() );
    }

    for( auto&& filename : filenames )
    {
      if( verbose )
        std::cerr << "* Processing " << filename << "...";

      SimplicialComplex K;
    }
  }

  // Output ------------------------------------------------------------
  //
  // After merging diagrams of corresponding dimensions, output all of
  // them in text format. Again, this is not the most efficient format
  // but it simplifies the remainder of the pipeline.

  for( auto&& filename : filenames )
  {
    auto&& diagrams = diagramCollection[ filename ];
    auto basename   = aleph::utilities::basename( filename );
    basename        = aleph::utilities::stem( basename );
    unsigned index  = 0;

    for( auto&& D : diagrams )
    {
      auto output = "/tmp/"
                  + basename
                  + "_d" + std::to_string( index++ )
                  + ".txt";

      std::ofstream out( output );
      out << D;
    }
  }
}
