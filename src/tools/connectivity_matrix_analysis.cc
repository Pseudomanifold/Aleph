/*
  This is a tool shipped by 'Aleph - A Library for Exploring Persistent
  Homology'.

  Its purpose is to analyse the persistent homology of connectivity
  matrices. The tool bears some semblance to the *network analysis*
  tools, but focuses specifically on data sets whose weights are an
  interpretable correlation measure.
*/

#include <aleph/persistenceDiagrams/Entropy.hh>
#include <aleph/persistenceDiagrams/Norms.hh>

#include <aleph/persistenceDiagrams/io/JSON.hh>

#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/io/AdjacencyMatrix.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/utilities/Filesystem.hh>

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <getopt.h>

void usage()
{
  std::cerr << "Usage: connectivity_matrix_analysis [--dimension DIMENSION] FILENAMES\n"
            << "\n"
            << "Analyses a set of connectivity matrices. The matrices are optionally\n"
            << "expanded to a pre-defined dimension. By default, only information of\n"
            << "the zeroth persistent homology group will be shown.\n"
            << "\n"
            << "Flags:\n"
            << "  -k: keep & report unpaired simplices (infinite values)\n"
            << "\n";
}

int main( int argc, char** argv )
{
  static option commandLineOptions[] =
  {
    { "dimension"     , required_argument, nullptr, 'd' },
    { "keep-unpaired" , no_argument      , nullptr, 'k' },
    { nullptr         , 0                , nullptr,  0  }
  };

  unsigned dimension = 0;
  bool keepUnpaired  = false;

  {
    int option = 0;

    while( ( option = getopt_long( argc, argv, "d:k", commandLineOptions, nullptr ) ) != -1 )
    {
      switch( option )
      {
      case 'd':
        dimension = static_cast<unsigned>( std::stoul(optarg) );
        break;
      case 'k':
        keepUnpaired = false;
        break;
      }
    }
  }

  if( (argc - optind) < 1 )
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

  for( int i = optind; i < argc; i++ )
    filenames.push_back( argv[i] );

  aleph::topology::io::AdjacencyMatrixReader reader;
  reader.setIgnoreNaNs();
  reader.setIgnoreZeroWeights();

  // Whew, is there *really* no better way of specifying this strategy
  // here as a qualified name?
  reader.setVertexWeightAssignmentStrategy(
      aleph::topology::io::AdjacencyMatrixReader::VertexWeightAssignmentStrategy::AssignZero
  );

  std::cout << "{\n"
            << "\"diagrams\": [\n";

  for( auto&& filename : filenames )
  {
    if( verbose )
      std::cerr << "* Processing " << filename << "...";

    SimplicialComplex K;
    reader( filename, K );

    K.sort();

    bool dualize                    = true;
    bool includeAllUnpairedCreators = keepUnpaired;

    auto diagrams
      = aleph::calculatePersistenceDiagrams( K,
                                             dualize,
                                             includeAllUnpairedCreators );

    if( verbose )
      std::cerr << "finished\n";

    auto basename
      = aleph::utilities::basename( filename );

    for( auto&& diagram : diagrams )
    {
      // Stores additional data about each persistence diagram in order
      // to make it easier to keep track of information.
      std::map<std::string, std::string> kvs;

      kvs["total_persistence_1"] = std::to_string( aleph::totalPersistence( diagram, 1.0 ) );
      kvs["total_persistence_2"] = std::to_string( aleph::totalPersistence( diagram, 2.0 ) );

      kvs["persistent_entropy"]  = std::to_string( aleph::persistentEntropy( diagram ) );

      aleph::io::writeJSON( std::cout, diagram, basename, kvs );
    }
  }

  std::cout << "\n"
            << "]\n"
            << "}\n";
}
