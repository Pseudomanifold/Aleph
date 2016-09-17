#include "io/Function.hh"
#include "representations/Vector.hh"

#include "Defaults.hh"

#include "BoundaryMatrix.hh"
#include "PersistencePairingCalculation.hh"

#include <iostream>
#include <string>

using DataType       = double;
using Index          = unsigned;
using Representation = aleph::representations::Vector<Index>;

using BoundaryMatrix = aleph::BoundaryMatrix<Representation>;

int main( int argc, char** argv )
{
  std::string filename;

  if( argc == 1 )
    return -1;

  if( argc >= 2 )
    filename = argv[1];

  BoundaryMatrix boundaryMatrix;
  std::vector<DataType> functionValues;

  aleph::io::loadFunction( filename,
                           boundaryMatrix,
                           functionValues );

  auto pairing =
    aleph::calculatePersistencePairing<aleph::defaults::ReductionAlgorithm>(
        boundaryMatrix );

  for( auto&& pair : pairing )
    std::cerr << pair.first << "--" << pair.second << "\n";
}
