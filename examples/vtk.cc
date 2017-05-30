/*
  This is an example file shipped by  'Aleph - A Library for Exporing
  Persistent Homology'.

  This example demonstrates how to load a VTK structured grid from
  a file, convert it into a simplicial complex, and calculate its
  persistent homology.

  Demonstrated classes:

    - aleph::PersistenceDiagram
    - aleph::topology::Simplex
    - aleph::topology::SimplicialComplex
    - aleph::topology::VTKStructuredGridReader

  Demonstrated functions:

    - aleph::calculatePersistenceDiagrams
    - aleph::utilities::extension
    - Simplicial complex sorting (for filtrations)

  Original author: Bastian Rieck
*/

#include "persistenceDiagrams/PersistenceDiagram.hh"

#include "topology/Simplex.hh"
#include "topology/SimplicialComplex.hh"

#include "topology/io/VTK.hh"

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


}
