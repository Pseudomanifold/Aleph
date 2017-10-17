/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  It permits the calculation of different entropy measures defined for
  persistence diagrams.

  Original author: Bastian Rieck
*/

#include <aleph/persistenceDiagrams/Entropy.hh>
#include <aleph/persistenceDiagrams/PersistenceDiagram.hh>

#include <aleph/persistenceDiagrams/io/Raw.hh>

#include <iostream>
#include <string>
#include <vector>

#include <cmath>

using DataType           = double;
using PersistenceDiagram = aleph::PersistenceDiagram<DataType>;

struct Input
{
  std::string filename;
  PersistenceDiagram persistenceDiagram;
};

void usage()
{
}

int main( int argc, char** argv )
{
  if( argc <= 1 )
  {
    usage();
    return -1;
  }

  std::vector<Input> inputs;
  inputs.reserve( argc );

  for( int i = 1; i < argc; i++ )
  {
    std::string filename = argv[i];

    std::cerr << "* Loading '" << filename << "'...";

    auto&& diagram = aleph::io::load<DataType>( filename );

    Input input = {
      filename,
      diagram
    };

    inputs.push_back( input );

    std::cerr << "finished\n";
  }

  for( auto&& input : inputs )
  {
    auto e_nn = aleph::nearestNeighbourAreaEntropy( input.persistenceDiagram );
    auto e_rg = aleph::gridEntropy( input.persistenceDiagram, 20 );

    std::cout << e_nn << "\t" << e_rg << "\n";
  }
}
