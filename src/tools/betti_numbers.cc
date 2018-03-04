/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  Its purpose is to calculate the zeroth and first Betti number of
  a data set given in the form of a graph in GML format. For every
  input filename, the program attempts to extract an ID. This will
  be used to assign Betti numbers in the output.

  This tool has been used to perform auxiliary calculations during
  molecular dynamics simulation, as described by:

  > Kai Sdeo, Bastian Rieck, Filip Sadlo
  > Visualization of Fullerene Fragmentation
  > IEEE Pacific Visualization Symposium 2018
*/

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/io/GML.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/utilities/Filesystem.hh>

#include <iostream>
#include <map>
#include <regex>
#include <string>
#include <vector>

int main( int argc, char** argv )
{
  if( argc <= 1 )
    return -1;

  bool verbose            = false;

  using DataType          = double;
  using VertexType        = unsigned short;
  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  std::vector<std::string> filenames;

  for( int i = 1; i < argc; i++ )
    filenames.push_back( argv[i] );

  aleph::topology::io::GMLReader reader;

  // Maps a data set ID to its corresponding Betti number. This is
  // required in order to generate a curve that measures how these
  // numbers change.
  std::map<unsigned, unsigned> id_to_betti;

  for( auto&& filename : filenames )
  {
    if( verbose )
      std::cerr << "* Processing " << filename << "...";

    SimplicialComplex K;
    reader( filename, K );

    K.sort();

    bool dualize                    = true;
    bool includeAllUnpairedCreators = true;

    auto diagrams
      = aleph::calculatePersistenceDiagrams( K,
                                             dualize,
                                             includeAllUnpairedCreators );

    if( verbose )
      std::cerr << "finished\n";

    auto basename
      = aleph::utilities::basename( filename );

    std::regex reDataSetSuffix( "\\D*([[:digit:]]+).*" );
    std::smatch matches;

    if( std::regex_match( basename, matches, reDataSetSuffix ) )
    {
      unsigned id = static_cast<unsigned>( std::stoull( matches[1] ) );

      if( diagrams.size() >= 2 )
        id_to_betti[id] = static_cast<unsigned>( diagrams[1].betti() );
      else
        id_to_betti[id] = static_cast<unsigned>( 0 );
    }
    else
      throw std::runtime_error( "Unable to identify ID" );
  }

  if( verbose )
    std::cerr << "* Obtained " << id_to_betti.size() << " data sets\n";

  for( auto&& pair : id_to_betti )
    std::cout << pair.first << "\t" << pair.second << "\n";

  std::cout << "\n\n";
}
