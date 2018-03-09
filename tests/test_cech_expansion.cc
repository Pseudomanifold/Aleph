#include <tests/Base.hh>

#include <aleph/geometry/CechComplex.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

#include <algorithm>
#include <iterator>
#include <set>
#include <vector>

#include <cmath>

using namespace aleph::geometry;
using namespace aleph::topology;
using namespace aleph;

template <class Data, class Vertex> void triangle()
{
  ALEPH_TEST_BEGIN( "Triangle" );

  ALEPH_TEST_END();
}

int main()
{
  triangle<double, unsigned>();
  triangle<double, short   >();
  triangle<float,  unsigned>();
  triangle<float,  short   >();
}
