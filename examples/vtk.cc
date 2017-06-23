/*
  This is an example file shipped by 'Aleph - A Library for Exploring
  Persistent Homology'.

  This example demonstrates how to load a VTK structured grid from
  a file, convert it into a simplicial complex, and calculate its
  persistent homology.

  Demonstrated classes:

    - aleph::PersistenceDiagram
    - aleph::topology::Simplex
    - aleph::topology::SimplicialComplex
    - aleph::topology::io::VTKStructuredGridReader

  Demonstrated functions:

    - aleph::calculatePersistenceDiagrams
    - aleph::utilities::extension
    - Simplicial complex sorting (for filtrations)

  Original author: Bastian Rieck
*/

#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/topology/io/VTK.hh>

#include <functional>
#include <iostream>
#include <string>

#include <getopt.h>

void usage()
{
  std::cerr << "Usage: vtk [--superlevels] [--sublevels] FILE\n"
            << "\n"
            << "Calculates persistent homology of an input file. This program only\n"
            << "handles VTK files. The output is a persistence diagram and will be\n"
            << "written to STDOUT.\n\n"
            << "Flags:\n"
            << "  -S --superlevels: calculate superlevel sets\n"
            << "  -s --sublevels  : calculate sublevel sets\n"
            << "\n";
}

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "superlevels", no_argument, nullptr, 'S' },
    { "sublevels"  , no_argument, nullptr, 's' },
    { nullptr      , 0          , nullptr,  0  }
  };

  bool calculateSuperlevelSets = false;

  int option = 0;
  while( ( option = getopt_long( argc, argv, "Ss", commandLineOptions, nullptr ) ) != -1 )
  {
    switch( option )
    {
    case 'S':
      calculateSuperlevelSets = true;
      break;
    case 's':
      calculateSuperlevelSets = false;
      break;
    default:
      break;
    }
  }

  if( ( argc - optind ) <= 0 )
  {
    usage();
    return -1;
  }

  std::string filename = argv[optind];

  // Since all data types in Aleph are templated, you need to decide
  // which types to use beforehand---or use more advanced constructs
  // like a type switch over a list of admissible types.
  using DataType           = double;
  using VertexType         = unsigned;
  using Simplex            = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex  = aleph::topology::SimplicialComplex<Simplex>;

  // Select a functor for calculating the weights when reading
  // a simplicial complex.
  //
  // This is an optional argument for the main operator of the
  // structured grid reader class. It is used to determine the
  // weight of an edge.
  //
  // Since sublevel sets 'grow' from small to large weights, a
  // correct assignment uses the `max()` function. Analogously
  // the `min()` function is used for superlevel sets.
  //
  // If you do not specify this functor, sublevel sets will be
  // calculated by default.
  auto functor = calculateSuperlevelSets
               ? [] ( DataType a, DataType b ) { return std::min(a,b); }
               : [] ( DataType a, DataType b ) { return std::max(a,b); };

  std::cerr << "* Loading '" << filename << "'...";

  SimplicialComplex K;
  aleph::topology::io::VTKStructuredGridReader reader;

  reader( filename, K, functor );

  std::cerr << "finished\n";

  std::cerr << "* Calculating persistent homology...";

  // Sort superlevel sets from largest to smallest weight...
  if( calculateSuperlevelSets )
    K.sort( aleph::topology::filtrations::Data<Simplex, std::greater<DataType> >() );
  // ...and sublevel sets from smallest to largest.
  else
    K.sort( aleph::topology::filtrations::Data<Simplex, std::less<DataType> >() );

  auto persistenceDiagrams = aleph::calculatePersistenceDiagrams( K );

  std::cerr << "finished\n";

  for( auto&& D : persistenceDiagrams )
  {
    // Remove all features with a persistence of zero. This is not
    // strictly required but it helps unclutter the diagram.
    D.removeDiagonal();

    std::cout << D << "\n";
  }
}
