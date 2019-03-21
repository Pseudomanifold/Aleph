/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  Its purpose is to analyse the persistent homology of connectivity
  matrices. The tool bears some semblance to the *network analysis*
  tools, but focuses specifically on data sets whose weights are an
  interpretable correlation measure.
*/

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/io/AdjacencyMatrix.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/utilities/Filesystem.hh>

#include <iostream>
#include <map>
#include <regex>
#include <string>
#include <vector>

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

  bool verbose            = false;

  using DataType          = double;
  using VertexType        = unsigned short;
  using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  std::vector<std::string> filenames;

  for( int i = 1; i < argc; i++ )
    filenames.push_back( argv[i] );

  aleph::topology::io::AdjacencyMatrixReader reader;

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

    // TODO: output of generated persistence diagrams, subject to
    // additional transformations
  }
}
