#include "config/Defaults.hh"

#include "persistentHomology/Calculation.hh"

#include "representations/Vector.hh"

#include "topology/BoundaryMatrix.hh"
#include "topology/io/Function.hh"

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

  auto diagram
    = aleph::calculatePersistenceDiagram<aleph::defaults::ReductionAlgorithm>(
        boundaryMatrix,
        functionValues );

  std::cerr << diagram << "\n";
}
