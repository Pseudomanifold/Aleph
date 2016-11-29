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
  RipsExpander<SimplicialComplex> ripsExpander;

  auto vr1 = ripsExpander( K, 2 );
  auto vr2 = ripsExpander( K, 3 );

  ALEPH_ASSERT_THROW( vr1.empty() == false );
  ALEPH_ASSERT_THROW( vr2.empty() == false );
  ALEPH_ASSERT_THROW( vr1.size()  == vr2.size() );
  ALEPH_ASSERT_THROW( vr1.size()  == 7 );

  ALEPH_TEST_END();
  return true;
}

int main()
{
  triangle<double, unsigned>();
}
