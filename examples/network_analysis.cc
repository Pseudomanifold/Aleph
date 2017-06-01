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
    - aleph::utilities::extension
    - Betti numbers of persistence diagrams
    - Simplicial complex sorting (for filtrations)

  Original author: Bastian Rieck
*/

#include "geometry/RipsExpander.hh"

#include "topology/io/EdgeLists.hh"
#include "topology/io/GML.hh"
#include "topology/io/Pajek.hh"

#include "persistentHomology/Calculation.hh"

#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include "topology/filtrations/Data.hh"

#include "utilities/Filesystem.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

void usage()
{
  std::cerr << "Usage: network_analysis FILE [DIMENSION]\n"
            << "\n"
            << "Loads a weighted network (graph) from FILE, expands it up to\n"
            << "the specified DIMENSION, and calculates persistence diagrams\n"
            << "of the weight function of the input.\n"
            << "\n"
            << "Diagrams will be written to STDOUT in a gnuplot-like style.\n"
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

  // These are just given for convenience. It makes declaring an
  // instance of a simplicial complex much easier.
  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  if( argc <= 1 )
  {
    usage();
    return -1;
  }

  // By convention, simplicial complexes are always called K (or L if
  // you need two of them). This is probably bad style. ;)
  SimplicialComplex K;

  // Reading -----------------------------------------------------------

  std::string filename = argv[1];
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

  aleph::geometry::RipsExpander<SimplicialComplex> ripsExpander;
  if( argc >= 3 )
    K = ripsExpander( K, static_cast<unsigned>( std::stoul( argv[2] ) ) );

  // This leaves the simplicial complex untouched. The simplices with
  // highest dimension are the edges, i.e. the 1-simplices.
  else
    K = ripsExpander( K, 1 );

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

    // This 'header' contains some informative entries about the
    // persistence diagram.
    //
    // The Betti number counts how many simplices are without a
    // partner in the diagram.
    //
    // Notice that the diagram has pre-defined output operator;
    // if you do not like the output, you may of course iterate
    // over the points in the diagram yourself.
    std::cout << "# Persistence diagram <" << filename << ">\n"
              << "#\n"
              << "# Dimension   : " << D.dimension() << "\n"
              << "# Entries     : " << D.size() << "\n"
              << "# Betti number: " << D.betti() << "\n"
              << D << "\n\n";
  }
}
