#include "complexes/FLANN.hh"
#include "complexes/VietorisRipsComplex.hh"

#include "containers/PointCloud.hh"

#include "distances/Euclidean.hh"

#include "utilities/String.hh"

#include <iostream>
#include <string>

using namespace aleph;
using namespace complexes;
using namespace geometry;
using namespace utilities;
using namespace distances;

int main( int argc, char** argv )
{
  if( argc <= 2 )
    return -1;

  using DataType   = double;
  using PointCloud = PointCloud<DataType>;
  using Distance   = Euclidean<DataType>;

  std::string input = argv[1];
  auto pointCloud   = load<DataType>( input );
  auto dimension    = pointCloud.dimension() + 1;
  auto epsilon      = convert<DataType>( argv[2] );

  if( argc >= 4 )
    dimension = std::stoul( argv[3] );

  std::cerr << "* Calculating Vietoris--Rips complex with eps=" << epsilon << " and d=" << dimension << "...";

  FLANN<PointCloud, Distance> flannWrapper( pointCloud );

  auto K
    = buildVietorisRipsComplex( flannWrapper,
                                epsilon,
                                unsigned( dimension ) );

  std::cerr << "finished\n"
            << "* Obtained simplicial complex with " << K.size() << " simplices\n";
}
