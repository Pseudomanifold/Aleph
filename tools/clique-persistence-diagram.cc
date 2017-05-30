#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>

#include <cassert>
#include <cmath>

// TODO: Replace this as soon as possible with a more modern option
// parser interface.
#include <getopt.h>

#include "geometry/RipsExpander.hh"
#include "geometry/RipsExpanderTopDown.hh"

#include "persistenceDiagrams/Norms.hh"
#include "persistenceDiagrams/PersistenceDiagram.hh"

#include "persistentHomology/ConnectedComponents.hh"

#include "topology/CliqueGraph.hh"
#include "topology/ConnectedComponents.hh"
#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include "topology/filtrations/Data.hh"

#include "topology/io/EdgeLists.hh"
#include "topology/io/GML.hh"
#include "topology/io/Pajek.hh"

#include "utilities/Filesystem.hh"

using DataType           = double;
using VertexType         = unsigned;
using Simplex            = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex  = aleph::topology::SimplicialComplex<Simplex>;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

namespace
{

// Basic functor for calculating all sorts of additional information
// about clique communities. This is easier than using the algorithm
// interface in order to provide said information.
class CliqueCommunityInformationFunctor
{
public:
  CliqueCommunityInformationFunctor( SimplicialComplex& K )
    : _K( K )
  {
    std::vector<VertexType> vertices;
    _K.vertices( std::inserter( vertices, vertices.end() ) );

    // Ensure proper initialization for all vertices, regardless of
    // whether they appear in any clique community or not.
    for( auto&& vertex : vertices )
    {
      _vim[vertex].accumulatedPersistence    = DataType();
      _vim[vertex].numberOfCliqueCommunities = 0;
    }
  }

  void initialize( VertexType v )
  {
    _cs[v] = 1;
    _cc[v] = {v};
  }

  void operator()( VertexType younger,  // younger connected component
                   VertexType older,    // older connected component
                   DataType creation,   // creation threshold
                   DataType destruction // destruction threshold
                 )
  {
    // Increase component size recursively by increasing the vertex
    // count for the older connected component.
    _cs[older] += _cs[younger];

    // Copy all vertices of the younger component to the older component
    // in order to signal that they are now merged.
    _cc[older].insert( _cc[older].end(),
                       _cc[younger].begin(), _cc[younger].end() );

    std::unordered_set<VertexType> cliqueVertices;

    for( auto&& vertex : _cc[younger] )
    {
      auto&& simplex = _K.at( vertex );

      cliqueVertices.insert( simplex.begin(),
                             simplex.end() );
    }

    // I don't want to count clique communities of negligible
    // persistence...
    if( creation != destruction )
    {
      for( auto&& vertex : cliqueVertices )
      {
        _vim[vertex].accumulatedPersistence    += std::pow( DataType( destruction - creation ), DataType(2) );
        _vim[vertex].numberOfCliqueCommunities += 1;
      }
    }

    _cc.erase(younger);
  }

  void operator()( VertexType root,   // root of an essential connected component
                   DataType creation  // creation threshold
                 )
  {
    std::unordered_set<VertexType> cliqueVertices;

    for( auto&& vertex : _cc[root] )
    {
      auto&& simplex = _K.at( vertex );

      cliqueVertices.insert( simplex.begin(),
                             simplex.end() );
    }

    for( auto&& vertex : cliqueVertices )
    {
      _vim[vertex].accumulatedPersistence    += std::pow( DataType( _destruction - creation ), DataType(2) );
      _vim[vertex].numberOfCliqueCommunities += 1;
    }
  }

  void setDestructionThreshold( DataType threshold )
  {
    _destruction = threshold;
  }

  /** Query component size information */
  unsigned getComponentSize( VertexType vertex ) const
  {
    return _cs.at( vertex );
  }

  /** Query accumulated persistence information */
  DataType accumulatedPersistence( VertexType vertex ) const
  {
    return _vim.at( vertex ).accumulatedPersistence;
  }

  /** Query number of clique communities */
  unsigned numberOfCliqueCommunities( VertexType vertex ) const
  {
    return _vim.at( vertex ).numberOfCliqueCommunities;
  }

private:

  // Original simplicial complex for looking up the vertices during
  // merging and centrality calculations.
  SimplicialComplex& _K;

  // Vertex information storage class. That way, I only require
  // a single map for looking up vertex information. The single
  // component map is not a mistake, though: it needs to remain
  // outside this struct because information needs to be erased
  // from it.
  struct VertexInformation
  {
    unsigned numberOfCliqueCommunities;
    DataType accumulatedPersistence;
  };

  std::unordered_map<VertexType, VertexInformation> _vim;

  // Maps with relative vertex indices --------------------------------
  //
  // The maps below use indices relative to the persistence diagram in
  // that dimension. Hence, they count the number of features or their
  // sizes without having direct knowledge about any other dimensions.

  std::unordered_map<VertexType, unsigned> _cs;                 // Component sizes
  std::unordered_map<VertexType, std::vector<VertexType> > _cc; // Connected components

  /**
    Destruction threshold to use for essential persistent homology
    classes. This threshold needs to be set by the client.
  */

  DataType _destruction = std::numeric_limits<DataType>::infinity();
};

} // anonymous namespace

std::string formatOutput( const std::string& prefix, unsigned k, unsigned K )
{
  std::ostringstream stream;
  stream << prefix;
  stream << std::setw( int( std::log10( K ) + 1 ) ) << std::setfill( '0' ) << k;
  stream << ".txt";

  return stream.str();
}

std::string formatLabel( const std::string label )
{
  // No whitespace---nothing to do
  if( label.find( '\t' ) == std::string::npos && label.find( ' ' ) == std::string::npos )
    return label;
  else
    return "\""+label+"\"";
}

void usage()
{
  // TODO: This is outdated...
  std::cerr << "Usage: clique-persistence-diagram [--invert-weights] [--reverse] FILE K\n"
            << "\n"
            << "Calculates the clique persistence diagram for FILE, which is\n"
            << "supposed to be a weighted graph. The K parameter denotes the\n"
            << "maximum dimension of a simplex for extracting a clique graph\n"
            << "and tracking persistence of clique communities.\n\n"
            << ""
            << "******************\n"
            << "Optional arguments\n"
            << "******************\n"
            << "\n"
            << " --centrality    : If specified, calculates centralities for\n"
            << "                   all vertices. Note that this uses copious\n"
            << "                   amounts of time because *all* communities\n"
            << "                   need to be extracted and inspected.\n"
            << "\n"
            << " --invert-weights: If specified, inverts input weights. This\n"
            << "                   is useful if the original weights measure\n"
            << "                   the strength of a relationship, and not a\n"
            << "                   dissimilarity.\n"
            << "\n"
            << " --reverse       : Reverses the enumeration order of cliques\n"
            << "                   by looking for higher-dimensional cliques\n"
            << "                   before enumerating lower-dimensional ones\n"
            << "                   instead of the other way around.\n"
            << "\n\n";
}

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "ignore-empty"  , no_argument      , nullptr, 'e' },
    { "invert-weights", no_argument      , nullptr, 'i' },
    { "normalize"     , no_argument      , nullptr, 'n' },
    { "reverse"       , no_argument      , nullptr, 'r' },
    { "min-k"         , required_argument, nullptr, 'k' },
    { nullptr         , 0                , nullptr,  0  }
  };

  bool ignoreEmpty         = false;
  bool invertWeights       = false;
  bool normalize           = false;
  bool reverse             = false;
  unsigned minK            = 0;

  int option = 0;
  while( ( option = getopt_long( argc, argv, "k:einr", commandLineOptions, nullptr ) ) != -1 )
  {
    switch( option )
    {
    case 'k':
      minK = unsigned( std::stoul( optarg ) );
      break;

    case 'e':
      ignoreEmpty = true;
      break;

    case 'i':
      invertWeights = true;
      break;

    case 'n':
      normalize = true;
      break;

    case 'r':
      reverse = true;
      break;

    default:
      break;
    }
  }

  if( (argc - optind ) < 2 )
  {
    usage();
    return -1;
  }

  std::string filename = argv[optind++];
  unsigned maxK        = static_cast<unsigned>( std::stoul( argv[optind++] ) );

  SimplicialComplex K;

  // Input -------------------------------------------------------------

  std::cerr << "* Reading '" << filename << "'...";

  // Optional map of node labels. If the graph contains node labels and
  // I am able to read them, this map will be filled.
  std::map<VertexType, std::string> labels;

  if( aleph::utilities::extension( filename ) == ".gml" )
  {
    aleph::topology::io::GMLReader reader;
    reader( filename, K );

    auto labelMap = reader.getNodeAttribute( "label" );

    // Note that this assumes that the labels are convertible to
    // numbers.
    //
    // TODO: Solve this generically?
    for( auto&& pair : labelMap )
      if( !pair.second.empty() )
        labels[ static_cast<VertexType>( std::stoul( pair.first ) ) ] = pair.second;

    if( labels.empty() )
      labels.clear();
  }
  else if( aleph::utilities::extension( filename ) == ".net" )
  {
    aleph::topology::io::PajekReader reader;
    reader( filename, K );

    auto labelMap = reader.getLabelMap();

    // Note that this assumes that the labels are convertible to
    // numbers.
    //
    // TODO: Solve this generically?
    for( auto&& pair : labelMap )
      if( !pair.second.empty() )
        labels[ static_cast<VertexType>( std::stoul( pair.first ) ) ] = pair.second;

    if( labels.empty() )
      labels.clear();

  }
  else
  {
    aleph::topology::io::EdgeListReader reader;
    reader.setReadWeights( true );
    reader.setTrimLines( true );

    reader( filename, K );
  }

  std::cerr << "finished\n";

  DataType maxWeight = std::numeric_limits<DataType>::lowest();
  DataType minWeight = std::numeric_limits<DataType>::max();
  for( auto&& simplex : K )
  {
    maxWeight = std::max( maxWeight, simplex.data() );
    minWeight = std::min( minWeight, simplex.data() );
  }

  if( normalize && maxWeight != minWeight )
  {
    std::cerr << "* Normalizing weights to [0,1]...";

    auto range = maxWeight - minWeight;

    for (auto it = K.begin(); it != K.end(); ++it )
    {
      if( it->dimension() == 0 )
        continue;

      auto s = *it;

      s.setData( ( s.data() - minWeight ) / range );
      K.replace( it, s );
    }

    maxWeight = DataType(1);
    minWeight = DataType(0);

    std::cerr << "finished\n";
  }

  if( invertWeights )
  {
    std::cerr << "* Inverting filtration weights...";

    for( auto it = K.begin(); it != K.end(); ++it )
    {
      if( it->dimension() == 0 )
        continue;

      auto s = *it;
      s.setData( maxWeight - s.data() );

      K.replace( it, s );
    }

    std::cerr << "finished\n";
  }

  // Expansion ---------------------------------------------------------

  std::cerr << "* Expanding simplicial complex to k=" << maxK << "...";

  if( reverse )
  {
    aleph::geometry::RipsExpanderTopDown<SimplicialComplex> ripsExpander;
    auto L = ripsExpander( K, maxK, minK );
    K      = ripsExpander.assignMaximumWeight( L, K );
  }
  else
  {
    aleph::geometry::RipsExpander<SimplicialComplex> ripsExpander;
    K = ripsExpander( K, maxK );
    K = ripsExpander.assignMaximumWeight( K );
  }

  std::cerr << "finished\n"
            << "* Expanded simplicial complex has " << K.size() << " simplices\n";

  K.sort( aleph::filtrations::Data<Simplex>() );

  // Stores the accumulated persistence of vertices. Persistence
  // accumulates if a vertex participates in a clique community.
  std::map<VertexType, double> accumulatedPersistenceMap;

  // Stores the number of clique communities a vertex is a part of.
  // I am using this only for debugging the algorithm.
  std::map<VertexType, unsigned> numberOfCliqueCommunities;

  std::vector<double> totalPersistenceValues;
  totalPersistenceValues.reserve( maxK );

  CliqueCommunityInformationFunctor ccif( K );
  ccif.setDestructionThreshold( 2 * maxWeight );

  // By traversing the clique graphs in descending order I can be sure
  // that a graph will be available. Otherwise, in case of a minimum k
  // parameter and a reverted expansion, only empty clique graphs will
  // be traversed.
  for( unsigned k = maxK; k >= 1; k-- )
  {
    std::cerr << "* Extracting " << k << "-cliques graph...";

    auto C
        = aleph::topology::getCliqueGraph( K, k );

    C.sort( aleph::filtrations::Data<Simplex>() );

    std::cerr << "finished\n";

    std::cerr << "* " << k << "-cliques graph has " << C.size() << " simplices\n";

    if( !ignoreEmpty && C.empty())
    {
      std::cerr << "* Stopping here because no further cliques for processing exist\n";
      break;
    }

    auto&& tuple = aleph::calculateZeroDimensionalPersistenceDiagram<Simplex, aleph::traits::PersistencePairingCalculation<aleph::PersistencePairing<VertexType> > >( C, ccif );
    auto&& pd    = std::get<0>( tuple );
    auto&& pp    = std::get<1>( tuple );

    pd.removeDiagonal();

    if( !C.empty() )
    {
      using namespace aleph::utilities;
      auto outputFilename = formatOutput( "/tmp/" + stem( basename( filename ) ) + "_k", k, maxK );

      std::cerr << "* Storing output in '" << outputFilename << "'...\n";

      std::transform( pd.begin(), pd.end(), pd.begin(),
                      [&maxWeight] ( const PersistenceDiagram::Point& p )
                      {
                        if( !std::isfinite( p.y() ) )
                          return PersistenceDiagram::Point( p.x(), 2 * maxWeight );
                        else
                          return PersistenceDiagram::Point( p );
                      } );

      std::ofstream out( outputFilename );
      out << "# Original filename: " << filename << "\n";
      out << "# k                : " << k        << "\n";

      {
        auto itPoint = pd.begin();
        for( auto itPair = pp.begin(); itPair != pp.end(); ++itPair, ++itPoint )
        {
          auto&& simplex = C.at( itPair->first );
          auto&& vertex  = *simplex.begin();

          out << itPoint->x() << "\t" << itPoint->y() << "\t" << ccif.getComponentSize( vertex ) << "\n";
        }
      }
    }
  }

  {
    using namespace aleph::utilities;
    auto outputFilename = "/tmp/" + stem( basename( filename ) ) + ".txt";

    std::cerr << "* Storing accumulated persistence values in '" << outputFilename << "'...\n";

    std::ofstream out( outputFilename );

    std::set<VertexType> vertices;
    K.vertices( std::inserter( vertices, vertices.end() ) );

    for( auto&& vertex : vertices )
    {
      out << vertex
          << "\t" << ccif.accumulatedPersistence( vertex )
          << "\t" << ccif.numberOfCliqueCommunities( vertex )
          << ( labels.empty() ? "" : "\t" + formatLabel( labels.at( vertex ) ) ) << "\n";
    }
  }
}
