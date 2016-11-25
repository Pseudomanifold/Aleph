#include "complexes/RipsExpander.hh"

#include "tests/Base.hh"

#include "Simplex.hh"
#include "SimplicialComplex.hh"

#include <vector>

using namespace aleph::complexes;
using namespace aleph;

template <class Data, class Vertex> bool triangle()
{
  ALEPH_TEST_BEGIN( "Triangle" );

  using Simplex           = Simplex<Data, Vertex>;
  using SimplicialComplex = SimplicialComplex<Simplex>;

  std::vector<Simplex> simplices
    = { {0}, {1}, {2}, {0,1}, {0,2}, {1,2} };

  SimplicialComplex K( simplices.begin(), simplices.end() );
  RipsExpander<Simplex> ripsExpander;

  auto vr = ripsExpander( K, 2 );

  ALEPH_ASSERT_THROW( vr.empty() == false );

  ALEPH_TEST_END();
  return true;
}

int main()
{
  triangle<double, unsigned>();
}
