#ifndef ALEPH_PERSISTENT_HOMOLOGY_CONNECTED_COMPONENTS_HH__
#define ALEPH_PERSISTENT_HOMOLOGY_CONNECTED_COMPONENTS_HH__

#include "persistenceDiagrams/PersistenceDiagram.hh"

#include "persistentHomology/PersistencePairing.hh"

#include "topology/SimplicialComplex.hh"
#include "topology/UnionFind.hh"

#include <algorithm>
#include <tuple>
#include <vector>

namespace aleph
{

/**
  Calculates zero-dimensional persistent homology, i.e. tracking of connected
  components, for a given simplicial complex. This is highly-efficient, as it
  only requires a suitable 'Union--Find' data structure.

  As usual, the function assumes that the simplicial complex is in filtration
  order, meaning that faces are preceded by their cofaces. The function won't
  check this, though!
*/

template <class Simplex>
  std::tuple<
    PersistenceDiagram<typename Simplex::DataType>,
    PersistencePairing<typename Simplex::VertexType>
  >
calculateZeroDimensionalPersistenceDiagram( const topology::SimplicialComplex<Simplex>& K )
{
  using DataType   = typename Simplex::DataType;
  using VertexType = typename Simplex::VertexType;

  using namespace topology;

  // Extract {0,1}-simplices -----------------------------------------
  //
  // Note that there is a range predicate for the simplicial complex class that
  // does essentially the same thing. However, the results of the query must be
  // consistent with respect to the filtration order.
  std::vector<Simplex> simplices;

  std::copy_if( K.begin(), K.end(),
                std::back_inserter( simplices ),
                [] ( const Simplex& s ) { return s.dimension() <= 1; } );

  SimplicialComplex<Simplex> S
    = SimplicialComplex<Simplex>( simplices.begin(), simplices.end() );

  std::vector<VertexType> vertices;

  for( auto&& s : simplices )
  {
    if( s.dimension() == 0 )
      vertices.push_back( *s.begin() );
  }

  UnionFind<VertexType> uf( vertices.begin(), vertices.end() );
  PersistenceDiagram<DataType> pd;
  PersistencePairing<VertexType> pp;

  for( auto&& simplex : S )
  {
    // Only edges can destroy a component
    if( simplex.dimension() == 1 )
    {
      // Prepare component destruction -------------------------------

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
      auto uIndex = S.index( Simplex( youngerComponent ) );
      auto vIndex = S.index( Simplex( olderComponent ) );

      // The younger component must have the _larger_ index as it is born
      // _later_ in the filtration.
      if( uIndex < vIndex )
      {
        std::swap( youngerComponent, olderComponent );
        std::swap( uIndex, vIndex );
      }

      auto creation    = S[uIndex].data();
      auto destruction = simplex.data();

      uf.merge( youngerComponent, olderComponent );

      pd.add( creation        , destruction    );
      pp.add( youngerComponent, olderComponent );
    }
  }

  // Store information about unpaired simplices ----------------------
  //
  // All components in the Union--Find data structure now correspond to
  // essential 0-dimensional homology classes of the input complex.

  std::vector<VertexType> roots;
  uf.roots( std::back_inserter( roots ) );

  for( auto&& root : roots )
  {
    auto creator = *S.find( Simplex( root ) );

    pd.add( creator.data() );
    pp.add( root           );
  }

  return std::make_tuple( pd, pp );
}

} // namespace aleph

#endif
