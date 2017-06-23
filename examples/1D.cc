#include <aleph/config/Defaults.hh>

#include <aleph/persistentHomology/Calculation.hh>


#include <aleph/topology/BoundaryMatrix.hh>
#include <aleph/topology/io/Function.hh>
#include <aleph/topology/representations/Vector.hh>

#include <iostream>
#include <string>

using DataType       = double;
using Index          = unsigned;
using Representation = aleph::topology::representations::Vector<Index>;
using BoundaryMatrix = aleph::topology::BoundaryMatrix<Representation>;

int main( int argc, char** argv )
{
  std::string filename;

  if( argc == 1 )
    return -1;

  if( argc >= 2 )
    filename = argv[1];

  BoundaryMatrix boundaryMatrix;
  std::vector<DataType> functionValues;

  aleph::topology::io::loadFunction( filename,
                                     boundaryMatrix,
                                     functionValues );

  auto diagram
    = aleph::calculatePersistenceDiagram<aleph::defaults::ReductionAlgorithm>(
        boundaryMatrix,
        functionValues );

  std::cerr << diagram << "\n";
}
