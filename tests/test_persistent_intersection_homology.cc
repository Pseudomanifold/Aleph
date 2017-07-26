#include <tests/Base.hh>

#include <aleph/persistentHomology/algorithms/Standard.hh>
#include <aleph/persistentHomology/PhiPersistence.hh>

#include <aleph/topology/Conversions.hh>
#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <unordered_map>
#include <vector>

template <class T> void test()
{
  ALEPH_TEST_BEGIN( "Persistent intersection homology: simple example" );

  using Simplex           = aleph::topology::Simplex<T>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  std::vector<Simplex> simplices =
  {
    {0},
    {1},
    {2},
    {3},
    {4},
    {0,1}, {0,3}, {0,4},
    {1,2}, {1,4},
    {2,3}, {2,4},
    {3,4},
    {0,3,4}, // A
    {1,2,4}, // B
    {2,3,4}, // C
    {0,1,4}  // E
  };

  std::map<Simplex, bool> phi;

  for( auto&& simplex : simplices )
  {
    if( simplex.contains(4) == false || simplex.dimension() == 2 )
      phi[simplex] = true;
    else
      phi[simplex] = false;
  }

  SimplicialComplex K( simplices.begin(), simplices.end() );
  SimplicialComplex L;

  std::size_t s    = 0;
  std::tie( L, s ) =
    aleph::partition( K,
                      [&phi] ( const Simplex& s )
                      {
                        return phi.at(s);
                      } );

  ALEPH_ASSERT_EQUAL( K.size(), L.size() );

  auto boundaryMatrix = aleph::topology::makeBoundaryMatrix( L, s );
  auto indexA         = L.index( {0,3,4} );

  using IndexType     = typename decltype(boundaryMatrix)::Index;
  auto columnA        = boundaryMatrix.getColumn( static_cast<IndexType>( indexA ) );

  ALEPH_ASSERT_EQUAL( columnA.size(), 3 );

  aleph::persistentHomology::algorithms::Standard algorithm;
  algorithm( boundaryMatrix );

  unsigned numAllowableChains    = 0;
  unsigned numAllowableTwoChains = 0;

  for( IndexType i = 0; i < boundaryMatrix.getNumColumns(); i++ )
  {
    auto lowestOne = boundaryMatrix.getMaximumIndex( i );
    if( lowestOne.second && lowestOne.first <= s )
    {
      ++numAllowableChains;

      if( L.at(i).dimension() == 2 )
        ++numAllowableTwoChains;
    }
  }

  ALEPH_ASSERT_THROW( numAllowableChains >= numAllowableTwoChains );
  ALEPH_ASSERT_EQUAL( numAllowableTwoChains, 1 );

  ALEPH_TEST_END();
}

template <class T> void testWedgeOfTwoCircles()
{
  ALEPH_TEST_BEGIN( "Persistent intersection homology: wedge of two circles" );

  using Simplex           = aleph::topology::Simplex<T>;
  using SimplicialComplex = aleph::topology::SimplicialComplex<Simplex>;

  SimplicialComplex K = {
    {0},
    {1},
    {2},
    {3},
    {4},
    {5},
    {6},
    {0,1}, {0,6},
    {1,2},
    {2,3}, {2,5}, {2,6},
    {3,4},
    {4,5}
  };

  SimplicialComplex X0 = { {2} };
  SimplicialComplex X1 = K;

  auto D1 = aleph::calculateIntersectionHomology( K, {X0,X1}, aleph::Perversity( {-1,-1} ) );
  auto D2 = aleph::calculateIntersectionHomology( K, {X0,X1}, aleph::Perversity( { 0, 0} ) );

  ALEPH_ASSERT_EQUAL( D1.size(), 1 );
  ALEPH_ASSERT_EQUAL( D2.size(), 2 );

  ALEPH_ASSERT_EQUAL( D1[0].betti(), 2 );

  // TODO: is this correct? In his Ph.D. thesis "Analyzing Stratified
  // Spaces Using Persistent Versions of Intersection and Local
  // Homology", Bendich states that this should be 0...
  ALEPH_ASSERT_EQUAL( D2[0].betti(), 1 );
  ALEPH_ASSERT_EQUAL( D2[1].betti(), 2 );

  ALEPH_TEST_END();
}

int main(int, char**)
{
  test<float> ();
  test<double>();

  testWedgeOfTwoCircles<float> ();
  testWedgeOfTwoCircles<double>();
}
