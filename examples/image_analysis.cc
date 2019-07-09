#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/topology/io/Matrix.hh>

#include <functional>
#include <iostream>
#include <string>

#include <getopt.h>

using DataType          = float;
using VertexType        = unsigned;
using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;
using Filtration        = aleph::topology::filtrations::Data<Simplex>;

int main( int argc, char** argv )
{
  bool useSuperlevelSets = false;

  {
    static option commandLineOptions[] = {
      { "sublevel"   , no_argument, nullptr, 's' },
      { "superlevel" , no_argument, nullptr, 'S' },
      { nullptr      , 0          , nullptr,  0  }
    };

    int option = -1;
    while( ( option = getopt_long( argc, argv, "sS", commandLineOptions, nullptr ) ) != - 1)
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
  }

  if( ( argc - optind ) < 1 )
    return -1;

  std::string filename = argv[1];

  SimplicialComplex K;

  aleph::topology::io::MatrixReader reader;
  reader( filename, K );

  if( useSuperlevelSets )
  {
    using Filtration = 
      aleph::topology::filtrations::Data<Simplex, std::greater<DataType> >;
  
    K.sort( Filtration() );
  }
  else
    K.sort( Filtration() );

  auto diagrams = aleph::calculatePersistenceDiagrams( K );

  for( auto&& D : diagrams )
  {
    D.removeDiagonal();

    std::cout << D << "\n";
  }
}
