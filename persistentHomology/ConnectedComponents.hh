#ifndef ALEPH_PERSISTENT_HOMOLOGY_CONNECTED_COMPONENTS_HH__
#define ALEPH_PERSISTENT_HOMOLOGY_CONNECTED_COMPONENTS_HH__

#include "persistenceDiagrams/PersistenceDiagram.hh"

#include "persistentHomology/PersistencePairing.hh"

#include "topology/SimplicialComplex.hh"
#include "topology/UnionFind.hh"

#include "utilities/EmptyFunctor.hh"

#include <algorithm>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace aleph
{

namespace traits
{

template <class Pairing> class PersistencePairingCalculation
{
public:

  PersistencePairingCalculation( Pairing& pairing )
    : _pairing( pairing )
  {
  }

  using IndexType = typename Pairing::IndexType;

  void add( IndexType u, IndexType v )
  {
    _pairing.add( u, v );
  }

  void add( IndexType u )
  {
    _pairing.add( u );
  }

private:
  Pairing& _pairing;
};

template <class Pairing> class NoPersistencePairingCalculation
{
public:

  NoPersistencePairingCalculation( Pairing& pairing )
    : _pairing( pairing )
  {
  }

  using IndexType = typename Pairing::IndexType;

  void add( IndexType /* u */, IndexType /* v */ ) const noexcept
  {
  }

  void add( IndexType /* u */ ) const noexcept
  {
  }

private:
  Pairing& _pairing;
};

class DiagonalElementCalculation
{
public:
  template <class T> bool operator()( T /* creation */, T /* destruction */ ) const noexcept
  {
    return true;
  }
};

class NoDiagonalElementCalculation
{
public:
  template <class T> bool operator()( T creation, T destruction ) const noexcept
  {
    return creation != destruction;
  }
};

} // namespace traits

/**
  Calculates zero-dimensional persistent homology, i.e. tracking of connected
  components, for a given simplicial complex. This is highly-efficient, as it
  only requires a suitable 'Union--Find' data structure.

  As usual, the function assumes that the simplicial complex is in filtration
  order, meaning that faces are preceded by their cofaces. The function won't
  check this, though!
*/

template <
  class Simplex,
  class PairingCalculationTraits = traits::NoPersistencePairingCalculation< PersistencePairing<typename Simplex::VertexType> >,
  class ElementCalculationTraits = traits::NoDiagonalElementCalculation,
  class Functor = aleph::utilities::EmptyFunctor
>
  std::tuple<
    PersistenceDiagram<typename Simplex::DataType>,
    PersistencePairing<typename Simplex::VertexType>,
    std::unordered_map<typename Simplex::VertexType, unsigned>
  >
calculateZeroDimensionalPersistenceDiagram( const topology::SimplicialComplex<Simplex>& K, Functor&& functor = Functor() )
{
  using DataType   = typename Simplex::DataType;
  using VertexType = typename Simplex::VertexType;

  using namespace topology;

  std::vector<VertexType> vertices;
  K.vertices( std::back_inserter( vertices ) );

  UnionFind<VertexType> uf( vertices.begin(), vertices.end() );
  PersistenceDiagram<DataType> pd;                               // Persistence diagram
  PersistencePairing<VertexType> pp;                             // Persistence pairing
  std::unordered_map<VertexType, unsigned> cs;                   // Component sizes
  std::unordered_map<VertexType, std::vector<VertexType> > cc;   // Connected components
  std::unordered_map<VertexType, DataType> ap;                   // Accumulated persistence

  PairingCalculationTraits ct( pp );
  ElementCalculationTraits et;

  for( auto&& vertex : vertices )
  {
    cs[vertex] = 1;
    cc[vertex] = { vertex };

    functor.initialize( vertex );
  }

  for( auto&& simplex : K )
  {
    // Only edges can destroy a component; we may safely skip any other
    // simplex with a different dimension.
    if( simplex.dimension() != 1 )
      continue;

    // Prepare component destruction -----------------------------------

    VertexType u = *( simplex.begin() );
    VertexType v = *( simplex.begin() + 1 );

    auto youngerComponent = uf.find( u );
    auto olderComponent   = uf.find( v );

    // If the component has already been merged by some other edge, we are
    // not interested in it any longer.
    if( youngerComponent == olderComponent )
      continue;

    // Ensures that the younger component is always the first component. A
    // component is younger if it its parent vertex precedes the other one
    // in the current filtration.
    auto uIndex = K.index( Simplex( youngerComponent ) );
    auto vIndex = K.index( Simplex( olderComponent ) );

    // The younger component must have the _larger_ index as it is born
    // _later_ in the filtration.
    if( uIndex < vIndex )
    {
      std::swap( youngerComponent, olderComponent );
      std::swap( uIndex, vIndex );
    }

    auto creation    = K[uIndex].data();
    auto destruction = simplex.data();

    uf.merge( youngerComponent, olderComponent );

    cs[olderComponent] += cs[youngerComponent];

    cc[olderComponent].insert( cc[olderComponent].end(),
                               cc[youngerComponent].begin(), cc[youngerComponent].end() );

    functor( youngerComponent,
             olderComponent,
             creation,
             destruction );

    for( auto&& vertex : cc[youngerComponent] )
      ap[vertex] += DataType( destruction - creation );

    cc.erase(youngerComponent);

    if( et( creation, destruction ) )
    {
      pd.add( creation                         , destruction                                   );
      ct.add( static_cast<VertexType>( uIndex ), static_cast<VertexType>( K.index( simplex ) ) );
    }
    else
      cs.erase( youngerComponent );
  }

  // Store information about unpaired simplices ------------------------
  //
  // All components in the Union--Find data structure now correspond to
  // essential 0-dimensional homology classes of the input complex.

  std::vector<VertexType> roots;
  uf.roots( std::back_inserter( roots ) );

  for( auto&& root : roots )
  {
    auto creator = *K.find( Simplex( root ) );

    pd.add( creator.data()                                );
    ct.add( static_cast<VertexType>( K.index( creator ) ) );

    functor( root,
             creator.data() );
  }

  return std::make_tuple( pd, pp, cs );
}

} // namespace aleph

#endif
