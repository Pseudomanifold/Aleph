/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  Its purpose is to calculate clique communities of a weighted network
  and convert them to JSON. Data will be written to STDOUT so that one
  can store it directly in a file. An optional threshold parameter can
  be used to filter cliques. This is useful when calculating auxiliary
  visualizations of a data set.

  This tool follows the publication:

    Clique Community Persistence: A Topological Visual Analysis Approach for Complex Networks
    Bastian Rieck, Ulderico Fugacci, Jonas Lukasczyk, Heike Leitte
    Submitted to IEEE Vis 2017

  Please see https://submanifold.github.io/Aleph/Rieck17d.html for more
  usage details.
*/

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <string>
#include <sstream>

#include <getopt.h>

#include <cmath>

#include <aleph/geometry/RipsExpander.hh>

#include <aleph/persistentHomology/ConnectedComponents.hh>

#include <aleph/topology/CliqueGraph.hh>
#include <aleph/topology/ConnectedComponents.hh>
#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/topology/io/EdgeLists.hh>
#include <aleph/topology/io/GML.hh>

#include <aleph/utilities/Filesystem.hh>

using DataType           = double;
using VertexType         = unsigned;
using Simplex            = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex  = aleph::topology::SimplicialComplex<Simplex>;

template <class Simplex> std::string formatSimplex( const Simplex& s, bool useLabels, const std::map<VertexType, std::string>& labels )
{
  std::ostringstream stream;
  stream << "[";

  for( auto it = s.begin(); it != s.end(); ++it )
  {
    if( it != s.begin() )
      stream << ",";

    if( useLabels )
      stream << "\"" << labels.at( *it ) << "\"";
    else
      stream << *it;
  }

  stream << "]";

  return stream.str();
}


void usage()
{
  std::cerr << "Usage: clique_communities_to_json FILE THRESHOLD K\n"
            << "\n"
            << "Extracts clique communities from FILE, which is supposed to be\n"
            << "a weighted graph. In the subsequent calculation, an edge whose\n"
            << "weight is larger than THRESHOLD will be ignored. K denotes the\n"
            << "maximum dimension of a simplex for the clique graph extraction\n"
            << "and the clique community calculation. This does not correspond\n"
            << "to the dimensionality of the clique. Hence, a parameter of K=2\n"
            << "will result in calculating 3-clique communities because all of\n"
            << "the 2-simplices have 3 vertices.\n"
            << "\n"
            << "Optional arguments:\n"
            << "\n"
            << " --label         : Use labels instead of indices to refer to\n"
            << "                   individual cliques. This is in particular\n"
            << "                   relevant for applications in which labels\n"
            << "                   are important, e.g. literature networks.\n"
            << "\n"
            << " --invert-weights: If specified, inverts input weights. This\n"
            << "                   is useful if the original weights measure\n"
            << "                   the strength of a relationship, and not a\n"
            << "                   dissimilarity.\n"
            << "\n"
            << " --normalize     : Normalizes all weights to [0,1]. Use this\n"
            << "                   to compare multiple networks.\n"
            << "\n\n";
}

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "labels"        , no_argument, nullptr, 'l' },
    { "normalize"     , no_argument, nullptr, 'n' },
    { "invert-weights", no_argument, nullptr, 'i' },
    { nullptr         , 0          , nullptr,  0  }
  };

  bool useLabels     = false;
  bool normalize     = false;
  bool invertWeights = false;

  int option = 0;
  while( ( option = getopt_long( argc, argv, "iln", commandLineOptions, nullptr ) ) != -1 )
  {
    switch( option )
    {
    case 'i':
      invertWeights = true;
      break;

    case 'l':
      useLabels = true;
      break;

    case 'n':
      normalize = true;
      break;

    default:
      break;
    }
  }

  if( (argc - optind) < 3 )
  {
    usage();
    return -1;
  }

  std::string filename = argv[optind++];
  double threshold     = std::stod( argv[optind++] );
  unsigned maxK        = static_cast<unsigned>( std::stoul( argv[optind++] ) );

  SimplicialComplex K;

  // Input -------------------------------------------------------------

  // Optional map of node labels. If the graph contains node labels and
  // I am able to read them, this map will be filled.
  std::map<VertexType, std::string> labels;

  std::cerr << "* Reading '" << filename << "'...";

  if( aleph::utilities::extension( filename ) == ".gml" )
  {
    aleph::topology::io::GMLReader reader;
    reader( filename, K );

    auto labelMap = reader.getNodeAttribute( "label" );

    // Note that this assumes that the labels are convertible to
    // numbers.
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

  // Determining weights -----------------------------------------------

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

  // Thresholding ------------------------------------------------------

  {
    std::cerr << "* Filtering input data to threshold epsilon=" << threshold << "...";

    std::vector<Simplex> simplices;

    std::remove_copy_if( K.begin(), K.end(), std::back_inserter( simplices ),
                         [&threshold] ( const Simplex& s )
                         {
                           return s.data() > threshold;
                         } );

    K = SimplicialComplex( simplices.begin(), simplices.end() );

    std::cerr << "finished\n";
  }

  // Expansion ---------------------------------------------------------

  aleph::geometry::RipsExpander<SimplicialComplex> ripsExpander;
  K = ripsExpander( K, maxK );
  K = ripsExpander.assignMaximumWeight( K );

  K.sort( aleph::topology::filtrations::Data<Simplex>() );

  std::cout << "{\n"
            << "  \"" << threshold << "\": {\n";

  for( unsigned k = 1; k <= maxK; k++ )
  {
    std::cerr << "* Extracting " << k << "-cliques graph...";

    auto C
        = aleph::topology::getCliqueGraph( K, k );

    C.sort( aleph::topology::filtrations::Data<Simplex>() );

    std::cerr << "finished\n";

    std::cerr << "* " << k << "-cliques graph has " << C.size() << " simplices\n";

    auto uf = aleph::topology::calculateConnectedComponents( C );

    std::set<VertexType> roots;
    uf.roots( std::inserter( roots, roots.begin() ) );

    std::cerr << "* " << k << "-cliques graph has " << roots.size() << " connected components\n";

    std::cout << "    \"" << (k+1) << "\": [\n";

    for( auto&& root : roots )
    {
      // The vertex IDs stored in the Union--Find data structure
      // correspond to the indices of the simplicial complex. It
      // thus suffices to map them back.
      std::set<VertexType> vertices;
      uf.get( root, std::inserter( vertices, vertices.begin() ) );

      std::vector<Simplex> simplices;

      std::transform( vertices.begin(), vertices.end(), std::back_inserter( simplices ),
                      [&K] ( VertexType v )
                      {
                        return K.at(v);
                      } );

      std::sort( simplices.begin(), simplices.end() );

      std::cout << "            [";

      for( auto it = simplices.begin(); it != simplices.end(); ++it )
      {
        if( it != simplices.begin() )
          std::cout << ",";

        std::cout << formatSimplex( *it, useLabels, labels );
      }

      std::cout << "]";

      if( root != *roots.rbegin() )
        std::cout << ",";

      std::cout << "\n";
    }

    std::cout << "    ]";

    if( k < maxK )
      std::cout << "  ,";

    std::cout << "\n";
  }

  std::cout << "  }\n"
            << "}\n";
}
