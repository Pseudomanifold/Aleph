#ifndef ALEPH_TOPOLOGY_FLOYD_WARSHALL_HH__
#define ALEPH_TOPOLOGY_FLOYD_WARSHALL_HH__

#include <limits>
#include <unordered_map>

#include <aleph/math/SymmetricMatrix.hh>

namespace aleph
{

namespace topology
{

/**
  Implements the Floyd--Warshall algorithm for a weighted simplicial
  complex. The algorithm calculates the matrix of pairwise distances
  between *all* nodes.

  @param K Simplicial complex
  @returns Matrix of distances. The indexing of the matrix follows
           the order in which the *vertices* of the simplicial are
           encountered.
*/

template <class SimplicialComplex> auto floydWarshall( const SimplicialComplex& K )
  -> aleph::math::SymmetricMatrix<
      typename SimplicialComplex::ValueType::DataType,
      typename SimplicialComplex::ValueType::VertexType>

{
  using Simplex    = typename SimplicialComplex::ValueType;
  using DataType   = typename Simplex::DataType;
  using VertexType = typename Simplex::VertexType;
  using Matrix     = aleph::math::SymmetricMatrix<DataType, VertexType>;

  // Set up vertex-to-index lookup table -------------------------------

  std::unordered_map<Simplex, VertexType> vertex_to_index;

  {
    VertexType index = VertexType();
    for( auto&& s : K )
    {
      if( s.dimension() == 0 )
        vertex_to_index[s] = index++;
    }
  }

  // Set up matrix -----------------------------------------------------
  //
  // First, all distances are initialized to either zero (self) or
  // infinity (all others). Next, edge weights of the complex will
  // be added to the matrix.

  auto n = vertex_to_index.size();

  Matrix M( static_cast<VertexType>( n ) );

  for( VertexType i = VertexType(0); i < VertexType(n); i++ )
  {
    M(i,i) = DataType(0);
    for( auto j = i+1; j < VertexType(n); j++ )
      M(i,j) = std::numeric_limits<DataType>::has_infinity ? std::numeric_limits<DataType>::infinity() : std::numeric_limits<DataType>::max();
  }

  for( auto&& s : K )
  {
    if( s.dimension() == 1 )
    {
      auto&& u  = s[0];
      auto&& v  = s[1];
      auto&& iu = vertex_to_index.at(u);
      auto&& iv = vertex_to_index.at(v);

      M(iu, iv) = s.data();
    }
  }

  for( VertexType k = VertexType(0); k < VertexType(n); k++ )
  {
    for( VertexType i = VertexType(0); i < VertexType(n); i++ )
    {
      for( auto j = i+1; j < VertexType(n); j++ )
      {
        if( M(i,j) > M(i,k) + M(k,j) )
          M(i,j) = M(i,k) + M(k,j);
      }
    }
  }

  return M;
}

} // namespace topology

} // namespace aleph

#endif
