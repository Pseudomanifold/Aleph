#ifndef ALEPH_TOPOLOGY_CONVERSION_HH__
#define ALEPH_TOPOLOGY_CONVERSION_HH__

#include <aleph/config/Defaults.hh>

#include <aleph/topology/BoundaryMatrix.hh>

#include <algorithm>
#include <unordered_map>

namespace aleph
{

namespace topology
{

/**
  Converts a simplicial complex into its boundary matrix representation.
  An optional index may be used to stop converting simplices whose index
  is larger than the specified maximum.

  If no maximum index is specified, the boundary matrices created by the
  function are suitable for (persistent) homology. If a maximum index is
  given, however, the matrices are particularly suitable for calculating
  (persistent) intersection homology.
*/

template <
  class Representation = aleph::defaults::Representation,
  class SimplicialComplex
> BoundaryMatrix<Representation> makeBoundaryMatrix( const SimplicialComplex& K, std::size_t max = 0 )
{
  using Simplex = typename SimplicialComplex::ValueType;
  using Index   = typename BoundaryMatrix<Representation>::Index;

  BoundaryMatrix<Representation> M;
  M.setNumColumns( static_cast<Index>( K.size() ) );

  // Prepare index map -------------------------------------------------
  //
  // The idea is to map simplices to their index within the filtration
  // in order to speed up the conversion process.

  std::unordered_map<Simplex, Index> simplex_to_index;

  {
    Index i = Index(0);

    for( auto&& simplex : K )
      simplex_to_index[simplex] = i++;
  }

  Index j = Index(0);

  for( auto&& itSimplex = K.begin(); itSimplex != K.end(); ++itSimplex )
  {
    std::vector<Index> column;
    column.reserve( itSimplex->size() );

    for( auto&& itBoundary = itSimplex->begin_boundary();
         itBoundary != itSimplex->end_boundary();
         ++itBoundary )
    {
      column.push_back( simplex_to_index[ *itBoundary ] );
    }

    M.setColumn( j, column.begin(), column.end() );

    ++j;

    if( max && j >= max )
      break;
  }

  return M;
}

} // namespace topology

} // namespace aleph

#endif
