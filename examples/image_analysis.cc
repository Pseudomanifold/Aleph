#include <aleph/persistentHomology/Calculation.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <aleph/topology/io/Matrix.hh>

#include <iostream>
#include <string>

using DataType          = float;
using VertexType        = unsigned;
using Simplex           = aleph::topology::Simplex<DataType, VertexType>;
using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;
using Filtration        = aleph::topology::filtrations::Data<Simplex>;

int main( int argc, char** argv )
{
  if( argc <= 1 )
    return -1;

  std::string filename = argv[1];

  SimplicialComplex K;

  aleph::topology::io::MatrixReader reader;
  reader( filename, K );

  K.sort( Filtration() );

  auto diagrams = aleph::calculatePersistenceDiagrams( K );

  for( auto&& D : diagrams )
  {
    D.removeDiagonal();

    std::cout << D << "\n";
  }
}
