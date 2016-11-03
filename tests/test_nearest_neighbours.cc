#include <complexes/FLANN.hh>
#include <complexes/NearestNeighbours.hh>

#include <containers/PointCloud.hh>

#include <vector>

#include <cassert>

using namespace aleph::complexes;
using namespace aleph;

int main()
{
  using Container = PointCloud<double>;
  Container container = load<double>( "Iris.txt" );

  assert( container.size() == 150 );
  assert( container.dimension() == 4 );
  assert( container.size() == 2);

  FLANN<Container> flannWrapper( container );
}
