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
    - Persistence diagram pruning

  Original author: Bastian Rieck
*/

#include "persistenceDiagrams/PersistenceDiagram.hh"
#include "persistenceDiagrams/Norms.hh"

#include "persistentHomology/Calculation.hh"

#include "topology/io/PLY.hh"

#include "utilities/Timer.hh"

#include <iostream>
#include <vector>

int main( int argc, char** argv )
{
  std::string filename;
  std::string property = "quality";

  if( argc == 1 )
  {
    // TODO: usage
    return -1;
  }

  if( argc >= 2 )
    filename = argv[1];

  if( argc >= 3 )
    property = argv[2];

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
  // TODO:
  //   - Expansion (higher-dimensional simplices)
  //   - Different filtrations (superlevel, sublevel)

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
    // Displays some statistics about the persistence diagrams. The
    // total persistence (with some power $p$) refers to the sum of
    // all persistence values, while the $p$-norm is merely a power
    // of the total persistence value.
    //
    // Both norms are useful for determining the amount of activity
    // within a data set.
    std::cerr << "* Total degree-1 persistence: " << aleph::totalPersistence( D, 1.0 ) << "\n"
              << "* Total degree-2 persistence: " << aleph::totalPersistence( D, 2.0 ) << "\n"
              << "* 1-norm:                     " << aleph::pNorm( D, 1.0 ) << "\n"
              << "* 2-norm:                     " << aleph::pNorm( D, 2.0 ) << "\n";
  }
}
