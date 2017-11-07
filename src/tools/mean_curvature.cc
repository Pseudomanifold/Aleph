#include <aleph/config/Eigen.hh>

#include <aleph/containers/PointCloud.hh>

#include <aleph/geometry/TangentSpace.hh>

#include <string>

using namespace aleph;
using namespace containers;
using namespace geometry;

using DataType   = double;
using PointCloud = PointCloud<DataType>;

int main( int argc, char** argv )
{
#if EIGEN_VERSION_AT_LEAST(3,3,0)
  if( argc <= 2 )
    return -1;

  std::string filename = argv[1];
  unsigned k           = unsigned( std::stoul( argv[2] ) );

  auto pc = load<DataType>( filename );

  TangentSpace ts;
  auto curvature = ts( pc, k );
  auto n         = pc.size();

  for( decltype(n) i = 0; i < n; i++ )
  {
    auto p = pc[i];
    auto c = curvature[i];

    for( auto&& x : p )
      std::cout << x << " ";

    std::cout << c << "\n";
  }
#endif
}
