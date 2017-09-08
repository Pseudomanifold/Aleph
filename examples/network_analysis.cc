/*
  This is an example file shipped by 'Aleph - A Library for Exploring
  Persistent Homology'.

  This example demonstrates how to load a network---a graph---from
  a variety of input files.

  The graph will subsequently be expanded to a simplicial complex and
  the persistent homology of the graph will be calculated.

  Demonstrated classes:

    - aleph::PersistenceDiagram
    - aleph::geometry::RipsExpander
    - aleph::topology::io::EdgeListReader
    - aleph::topology::io::GMLReader
    - aleph::topology::io::PajekReader
    - aleph::topology::Simplex
    - aleph::topology::SimplicialComplex
    - aleph::topology::filtrations::Data

  Demonstrated functions:

    - aleph::calculatePersistenceDiagrams
    - aleph::utilities::basename
    - aleph::utilities::extension
    - aleph::utilities::format
    - aleph::utilities::stem
    - Betti numbers of persistence diagrams
    - Simplicial complex dimensionality queries
    - Simplicial complex sorting (for filtrations)

  Original author: Bastian Rieck
*/

#include <aleph/geometry/RipsExpander.hh>

#include <aleph/topology/io/EdgeLists.hh>
#include <aleph/topology/io/GML.hh>
#include <aleph/topology/io/Pajek.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/utilities/Filesystem.hh>
#include <aleph/utilities/Format.hh>

#include <cmath>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <getopt.h>

void usage()
{
  std::cerr << "Usage: network_analysis FILE [DIMENSION]\n"
            << "\n"
            << "Loads a weighted network (graph) from FILE, expands it up to\n"
            << "the specified DIMENSION, and calculates persistence diagrams\n"
            << "of the weight function of the input.\n"
            << "\n"
            << "Diagrams will be written to STDOUT in a gnuplot-like style.\n"
            << "\n"
            << "Optional arguments:\n"
            << "\n"
            << " --infinity FACTOR: Sets the value to use for unpaired points\n"
            << "                   in the persistence diagram. By default, a\n"
            << "                   large number or +inf will be used. If the\n"
            << "                   specified number is non-zero, it shall be\n"
            << "                   used as a factor in the weight assignment\n"
            << "                   of these points.\n"
            << "\n"
            << " --invert-weights: If specified, inverts input weights. This\n"
            << "                   is useful if the original weights measure\n"
            << "                   the strength of a relationship, and not a\n"
            << "                   dissimilarity between nodes.\n"
            << "\n"
            << " --node-weights  : Specifies a file from which to load node\n"
            << "                   weights for the filtration.\n"
            << "\n"
            << " --normalize     : Normalizes all weights to [0,1]. Use this\n"
            << "                   to compare multiple networks.\n"
            << "\n"
            << " --output PATH   : Uses the specified path to store diagrams\n"
            << "                   instead of writing them to STDOUT.\n"
            << "\n"
            << "\n";
}

int main( int argc, char** argv )
{
  // We have to specify the required data type and vertex type for the
  // simplex class here. These settings are generic and should cover
  // most situations. If your weights are floats or even integers, you
  // can of course change it here. Likewise, if your simplicial complex
  // is very small, you could use `short` instead of `unsigned` to
  // represent all possible values for vertices.
  using DataType          = double;
  using VertexType        = unsigned;

  static option commandLineOptions[] =
  {
    { "infinity"      , required_argument, nullptr, 'f' },
    { "invert-weights", no_argument      , nullptr, 'i' },
    { "node-weights"  , required_argument, nullptr, 'w' },
    { "normalize"     , no_argument      , nullptr, 'n' },
    { "output"        , required_argument, nullptr, 'o' },
    { nullptr         , 0                , nullptr,  0  }
  };

  bool invertWeights          = false;
  bool normalize              = false;
  DataType infinity           = std::numeric_limits<DataType>::has_infinity ? std::numeric_limits<DataType>::infinity() : std::numeric_limits<DataType>::max();
  std::string basePath        = std::string();
  std::string nodeWeightsFile = std::string();

  {
    int option = 0;
    while( ( option = getopt_long( argc, argv, "f:iw:no:", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'f':
        infinity = static_cast<DataType>( std::stod( optarg ) );
        break;

      case 'i':
        invertWeights = true;
        break;

      case 'n':
        normalize = true;
        break;

      case 'o':
        basePath = optarg;
        break;

      case 'w':
        nodeWeightsFile = optarg;
        break;

      default:
        break;
      }
    }
  }

  if( (argc - optind ) < 1 )
  {
    usage();
    return -1;
  }

  // These are just given for convenience. It makes declaring an
  // instance of a simplicial complex much easier.
  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  // By convention, simplicial complexes are always called K (or L if
  // you need two of them). This is probably bad style. ;)
  SimplicialComplex K;

  // Reading -----------------------------------------------------------

  std::string filename = argv[optind++];
  std::cerr << "* Reading '" << filename << "'...";

  // For supporting different graph formats, Aleph offers different
  // loader classes. Every loader class uses the same interface for
  // reading simplicial complexes, via `operator()()`.
  //
  // Some loader classes support special features of a given format
  // such as label reading. Here, we do not make use of any of them
  // because we are only interested in demonstrating the expansion.
  //
  // Every reader attempts to read edge weights from the file. They
  // are assigned to the 1-simplices in the simplicial complex.
  {
    auto extension = aleph::utilities::extension( filename );
    if( extension == ".gml" )
    {
      aleph::topology::io::GMLReader reader;
      reader( filename, K );
    }
    else if( extension == ".net" )
    {
      aleph::topology::io::PajekReader reader;
      reader( filename, K );
    }
    else
    {
      aleph::topology::io::EdgeListReader reader;
      reader( filename, K );
    }
  }

  std::cerr << "finished\n"
            << "* Extracted simplicial complex with " << K.size() << " simplices\n";

  // Assign nodes weights if specified by the user. This requires
  // re-calculating *all* weights of the simplicial complex.

  if( !nodeWeightsFile.empty() )
  {
    std::ifstream in( nodeWeightsFile );

    std::vector<DataType> nodeWeights;

    aleph::geometry::RipsExpander<SimplicialComplex> ripsExpander;
    K = ripsExpander.assignMaximumData(
      K,
      std::istream_iterator<DataType>( in ), std::istream_iterator<DataType>() );
  }

  // Pre-processing ----------------------------------------------------
  //
  // Determine the minimum and the maximum weight. If desired by the
  // user, normalize those weights and/or invert them.

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

  // Rips expansion ----------------------------------------------------
  //
  // The Rips expander is a generic class for expanding a graph into a
  // higher-dimensional simplicial complex. This follows the procedure
  // described in:
  //
  //   Fast construction of the Vietoris--Rips complex
  //   Afra Zomorodian
  //   Computers & Graphics, Volume 34, Issue 3, pp. 263-–271
  //
  // Formally speaking, the algorithm creates a k-simplex for subsets
  // of simplices in which k+1 simplices have pairwise intersections.
  //
  // Hence, if all three edges of a triangle are present, a 2-simplex
  // will be added to the simplicial complex.

  std::cerr << "* Expanding simplicial complex...";

  // By default, tThis leaves the simplicial complex untouched. The
  // simplices with highest dimension are the 1-simplices, i.e. the
  // edges. If the user specified an optional parameter, we use it.
  unsigned k = 1;

  if( argc - optind > 0 )
    k = static_cast<unsigned>( std::stoul( argv[optind++] ) );

  aleph::geometry::RipsExpander<SimplicialComplex> ripsExpander;
  K = ripsExpander( K, k );

  // This tells the expander to use the maximum weight of the faces of
  // a simplex in order to assign the weight of the simplex. Thus, the
  // simplicial complex models a sublevel set filtration after sorting
  // it with the appropriate functor.
  K = ripsExpander.assignMaximumWeight( K );

  std::cerr << "...finished\n"
            << "* Expanded complex has dimension " << K.dimension() << "\n"
            << "* Expanded complex has " << K.size() << " simplices\n";

  std::cerr << "* Establishing filtration order...";

  // A filtration is an ordering of K in which faces precede co-faces in
  // order to ensure that the ordering represents a growing complex. The
  // functor used below sorts simplices based on their data. It hence is
  // just another way of expressing the weights specified in the graph.
  //
  // The data-based filtration ensures that upon coinciding weights, the
  // lower-dimensional simplex will precede the higher-dimensional one.
  K.sort(
    aleph::topology::filtrations::Data<Simplex>()
  );

  std::cerr << "...finished\n";

  // Persistent homology calculation -----------------------------------
  //
  // This uses a convenience function. It expects a simplicial complex
  // in filtration order, calculates it persistent homology using the
  // default algorithms, and returns a vector of persistence diagrams.

  std::cerr << "* Calculating persistent homology...";

  auto diagrams
    = aleph::calculatePersistenceDiagrams( K );

  std::cerr << "...finished\n";

  for( auto&& D : diagrams )
  {
    // This removes all features with a persistence of zero. They only
    // clutter up the diagram.
    D.removeDiagonal();

    if( std::isfinite( infinity ) && infinity != std::numeric_limits<DataType>::max() )
    {
      std::cerr << "* Transforming unpaired points in persistence diagram with a factor of " << infinity << "...\n";

      using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

      std::transform( D.begin(), D.end(), D.begin(),
          [&maxWeight, &infinity] ( const PersistenceDiagram::Point& p )
          {
            if( !std::isfinite( p.y() ) )
              return PersistenceDiagram::Point( p.x(), infinity * maxWeight );
            else
              return PersistenceDiagram::Point( p );
          } );
    }

    std::ostringstream stream;

    // This 'header' contains some informative entries about the
    // persistence diagram.
    //
    // The Betti number counts how many simplices are without a
    // partner in the diagram.
    //
    // Notice that the diagram has pre-defined output operator;
    // if you do not like the output, you may of course iterate
    // over the points in the diagram yourself.
    stream << "# Persistence diagram <" << filename << ">\n"
           << "#\n"
           << "# Dimension   : " << D.dimension() << "\n"
           << "# Entries     : " << D.size() << "\n"
           << "# Betti number: " << D.betti() << "\n"
           << D;

    // Use the separator for all but the last persistence diagram. Else,
    // the output format is inconsistent and results in empty files.
    if( D != diagrams.back() )
      stream << "\n\n";

    if( basePath.empty() )
      std::cout << stream.str();
    else
    {
      using namespace aleph::utilities;

      // The formatting ensures that if the user specified an expansion
      // up to, say, 10, the resulting files will use two-digit numbers
      // instead of one-digit numbers, regardless of how many complexes
      // there are.
      auto outputFilename = basePath + "/" + stem( basename( filename ) )
                                     + "_d"
                                     + format( D.dimension(), std::max( K.dimension(), static_cast<std::size_t>(k) ) )
                                     + ".txt";

      std::cerr << "* Storing output in '" << outputFilename << "'...\n";

      std::ofstream out( outputFilename );
      out << stream.str();
    }
  }
}
