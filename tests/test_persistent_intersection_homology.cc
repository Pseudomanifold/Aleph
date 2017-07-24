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

int main(int, char**)
{
  test<float> ();
  test<double>();
}
