/*
  This is an example file shipped by 'Aleph - A Library for Exploring
  Persistent Homology'.

  This example demonstrates how to load a mesh in PLY format from
  a file, convert it into a simplicial complex, and calculate its
  persistent homology. The Weights for the simplicial complex are
  taken from a user-specified property within the PLY file.

  Demonstrated classes:

    - aleph::PersistenceDiagram
    - aleph::topology::Simplex
    - aleph::topology::SimplicialComplex
    - aleph::topology::io::PLYReader
    - aleph::utilities::Timer

  Demonstrated functions:

    - aleph::calculatePersistenceDiagrams
    - aleph::pNorm
    - aleph::totalPersistence
    - Filtrations of simplicial complexes
    - Persistence diagram pruning (diagonal & unpaired values)

  Original author: Bastian Rieck
*/

#include "persistenceDiagrams/PersistenceDiagram.hh"
#include "persistenceDiagrams/Norms.hh"

#include "persistentHomology/Calculation.hh"

#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include "topology/filtrations/Data.hh"
#include "topology/io/PLY.hh"

#include "utilities/Timer.hh"

#include <iostream>
#include <vector>

#include <getopt.h>

void usage()
{
  std::cerr << "Usage: ply [--sublevel | --superlevel] FILENAME [PROPERTY]\n"
            << "\n"
            << "Reads a PLY mesh from FILENAME and converts it into a simplicial\n"
            << "complex. If specified, reads PROPERTY for each vertex (a quality\n"
            << "value, for example), and uses it to assign simplex weights.\n"
            << "\n"
            << "By default, the reader just uses the z coordinate of vertices in\n"
            << "the mesh because this property is guaranteed to exist.\n"
            << "\n"
            << "You may select a filtration for persistent homology calculations\n"
            << "using '--sublevel' (default) or '--superlevel'. This will change\n"
            << "the ordering of the simplicial complex, and thus the persistence\n"
            << "diagram.\n"
            << "\n"
            << "Flags:\n"
            << "  -s: use sublevel set filtration (default)\n"
            << "  -S: use superlevel set filtration\n"
            << "\n";
}

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "sublevel"   , no_argument, nullptr, 's' },
    { "superlevel" , no_argument, nullptr, 'S' },
    { nullptr      , 0          , nullptr,  0  }
  };

  bool useSuperlevelSets = false;

  int option = 0;
  while( ( option = getopt_long( argc, argv, "sS", commandLineOptions, nullptr ) ) != -1 )
  {
    switch( option )
    {
    case 's':
      useSuperlevelSets = false;
      break;
    case 'S':
      useSuperlevelSets = true;
      break;
    default:
      break;
    }
  }

  std::string filename;
  std::string property;

  if( ( argc - optind ) < 1 )
  {
    usage();
    return -1;
  }

  if( optind < argc )
    filename = argv[optind++];

  if( optind < argc )
    property = argv[optind++];

  // Loading -----------------------------------------------------------

  // We first declare the data types we want to convert the file to so
  // that the PLY reader class knows the desired type of complex. Note
  // that the `Simplex` and `SimplicialComplex` type are declared just
  // for reasons of convenience.
  using DataType          = double;
  using VertexType        = unsigned;
  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  // Declares the reader for loading a PLY file. At present, the loading
  // of binary files is not yet supported (even though some code exists)
  // and the reader only handles ASCII files properly.
  //
  // Note that we set a 'data property'. This specifies the attribute of
  // every vertex that is used to assign the data values of simplices in
  // the simplicial complex. The property defaults to 'z', i.e. the last
  // coordinate of every vertex. In this example, we permit users to use
  // another property.
  //
  // See https://en.wikipedia.org/wiki/PLY_(file_format) for information
  // about the PLY format.
  aleph::topology::io::PLYReader plyReader;

  if( !property.empty() )
    plyReader.setDataProperty( property );

  SimplicialComplex K;
  plyReader( filename, K );

  std::cerr << "* Loaded simplicial complex with " << K.size() << " simplices\n";

  // Persistent homology -----------------------------------------------
  //
  // At present, only the mesh connectivity is used to calculate
  // persistent homology.
  //
  // Note that we do not have to sort K, the simplicial complex,
  // prior to the calculations. By default, the PLY reader sorts
  // the complex from small weights to large weights.
  //
  // If the superlevel set filtration is requested, however, all
  // weights need to be re-calculated and the proper order needs
  // to be restored.
  if( useSuperlevelSets )
  {
    K.recalculateWeights( false );
    K.sort( aleph::topology::filtrations::Data<Simplex>() );
  }

  // This small utility class permits measuring the time of certain
  // operations. Internally, it makes use of `std::chrono` in order
  // to permit a sufficiently fine resolution.
  aleph::utilities::Timer timer;

  auto diagrams
    = aleph::calculatePersistenceDiagrams( K );

  std::cerr << "* Calculated " << diagrams.size() << " persistence diagrams in " << timer.elapsed_s() << "s\n";

  for( auto&& D : diagrams )
  {
    // Removes all features with zero persistence from the diagram in
    // order to simplify it.
    D.removeDiagonal();

    // This results in a tabular output of all points in the diagram;
    // it is ideally suited for further analysis in 'gnuplot'. Notice
    // that the two new lines can be used to automatically detect the
    // next dimension of the diagram.
    std::cout << D << "\n\n";
  }

  for( auto&& D : diagrams )
  {
    // Remove unpaired simplices of every persistence diagram prior
    // to calculating its statistics. Else, the norm will always be
    // infinite.
    D.removeUnpaired();

    // Displays some statistics about the persistence diagrams. The
    // total persistence (with some power $p$) refers to the sum of
    // all persistence values, while the $p$-norm is merely a power
    // of the total persistence value.
    //
    // Both norms are useful for determining the amount of activity
    // within a data set.
    std::cerr << "Dimension [" << D.dimension() << "]\n"
              << "* Total degree-1 persistence: " << aleph::totalPersistence( D, 1.0 ) << "\n"
              << "* Total degree-2 persistence: " << aleph::totalPersistence( D, 2.0 ) << "\n"
              << "* 1-norm:                     " << aleph::pNorm( D, 1.0 ) << "\n"
              << "* 2-norm:                     " << aleph::pNorm( D, 2.0 ) << "\n"
              << "\n";
  }
}
