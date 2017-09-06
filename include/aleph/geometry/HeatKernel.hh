#ifndef ALEPH_GEOMETRY_HEAT_KERNEL_HH__
#define ALEPH_GEOMETRY_HEAT_KERNEL_HH__

#include <aleph/config/Eigen.hh>

#ifdef ALEPH_WITH_EIGEN
  #include <Eigen/Core>
#endif

#include <unordered_map>
#include <vector>

namespace aleph
{

namespace geometry
{

#ifdef ALEPH_WITH_EIGEN

/**
  Extracts a weighted adjacency matrix from a simplicial complex. At
  present, this function only supports adjacencies between edges, so
  the resulting matrix is a graph adjacency matrix.

  @param K Simplicial complex

  @returns Weighted adjacency matrix. The indices of rows and columns
           follow the order of the vertices in the complex.
*/

template <class SimplicialComplex> auto weightedAdjacencyMatrix( const SimplicialComplex& K ) -> Eigen::Matrix<typename SimplicialComplex::ValueType::DataType, Eigen::Dynamic, Eigen::Dynamic>
{
  using Simplex    = typename SimplicialComplex::ValueType;
  using VertexType = typename Simplex::VertexType;
  using DataType   = typename Simplex::DataType;
  using Matrix     = Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>;

#if EIGEN_VERSION_AT_LEAST(3,3,0)
  using IndexType  = Eigen::Index;
#else
  using IndexType  = typename Matrix::Index;
#endif

  // Prepare map from vertex to index ----------------------------------

  std::unordered_map<VertexType, IndexType> vertex_to_index;
  IndexType n = IndexType();

  {
    std::vector<VertexType> vertices;
    K.vertices( std::back_inserter( vertices ) );

    IndexType index = IndexType();

    for( auto&& vertex : vertices )
      vertex_to_index[vertex] = index++;

    n = static_cast<IndexType>( vertices.size() );
  }

  // Prepare matrix ----------------------------------------------------

  Matrix W( n, n );

  for(auto&& s : K )
  {
    if( s.dimension() != 1 )
      continue;

    auto&& u = s[0];
    auto&& v = s[1];
    auto&& i = vertex_to_index.at( u );
    auto&& j = vertex_to_index.at( v );

    W(i,j)   = s.data();
    W(j,i)   = W(i,j);
  }

  return W;
}

template <class SimplicialComplex> auto weightedLaplacianMatrix( const SimplicialComplex& K ) -> Eigen::Matrix<typename SimplicialComplex::ValueType::DataType, Eigen::Dynamic, Eigen::Dynamic>
{
  auto W          = weightedAdjacencyMatrix( K );
  using Matrix    = decltype(W);
  using IndexType = typename Matrix::Index;

  Matrix L( W.rows(), W.cols() );

  auto V = W.rowwise().sum();

  for( IndexType i = 0; i < V.size(); i++ )
    L(i,i) = V(i);

  return L - W;
}

#endif

} // namespace geometry

} // namespace aleph

#endif
