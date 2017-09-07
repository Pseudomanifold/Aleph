#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/BetaSkeleton.hh>
#include <aleph/geometry/HeatKernel.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>

using DataType   = double;
using Distance   = aleph::distances::Euclidean<DataType>;
using PointCloud = aleph::containers::PointCloud<DataType>;

int main( int argc, char** argv )
{
  if( argc <= 1 )
    return -1;

  auto filename   = std::string( argv[1] );
  auto pointCloud = aleph::containers::load<DataType>( filename );

  std::cerr << "* Loaded point cloud with " << pointCloud.size() << " points\n";

  // Skeleton construction ---------------------------------------------

  // TODO: make configurable
  DataType beta = 1.0;

  std::cerr << "* Calculating beta-skeleton with beta = " << beta << "...";

  auto betaSkeleton
    = aleph::geometry::buildBetaSkeletonNaive( pointCloud,
                                               beta,
                                               Distance() );

  std::cerr << "...finished\n"
            << "* Simplical complex has " << betaSkeleton.size() << " simplices\n";

  // Scale estimation --------------------------------------------------
}
