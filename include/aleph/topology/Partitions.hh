#ifndef ALEPH_TOPOLOGY_PARTITIONS_HH__
#define ALEPH_TOPOLOGY_PARTITIONS_HH__

#include <aleph/config/Eigen.hh>

#ifdef ALEPH_WITH_EIGEN
  #include <Eigen/Core>
  #include <Eigen/Eigenvalues>
#endif

#include <aleph/geometry/HeatKernel.hh>

#include <aleph/math/Quantiles.hh>

#include <unordered_map>
#include <stdexcept>
#include <vector>

namespace aleph
{

namespace topology
{

template <class SimplicialComplex> std::vector<SimplicialComplex> bisect( const SimplicialComplex& K )
{
#ifdef ALEPH_WITH_EIGEN

  auto L = aleph::geometry::weightedLaplacianMatrix( K );

  Eigen::SelfAdjointEigenSolver< decltype(L) > solver;
  solver.compute( L );

  using Simplex    = typename SimplicialComplex::ValueType;
  using VertexType = typename Simplex::VertexType;
  using DataType   = typename Simplex::DataType;

  auto&& eigenvectors = solver.eigenvectors().template cast<DataType>();

  if( eigenvectors.size() < 2 )
    throw std::runtime_error( "Laplacian matrix dimensions are insufficient for bisection" );

  std::vector<DataType> fiedlerVector;

  {
    auto fiedlerVector_ = eigenvectors.col(1);

    fiedlerVector.assign( fiedlerVector_.data(),
                          fiedlerVector_.data() + fiedlerVector_.size() );
  }

  auto median     = aleph::math::median( fiedlerVector.begin(), fiedlerVector.end() );
  using IndexType = typename std::vector<DataType>::size_type;

  // Prepare map from index to vertex ----------------------------------

  std::unordered_map<IndexType, VertexType> index_to_vertex;

  {
    std::vector<VertexType> vertices;
    K.vertices( std::back_inserter( vertices ) );

    IndexType index = IndexType();

    for( auto&& vertex : vertices )
      index_to_vertex[index++] = vertex;
  }

  // Partition vertices ------------------------------------------------

  std::unordered_map<VertexType, bool> partition;

  for( IndexType i = 0; i < fiedlerVector.size(); i++ )
  {
    auto vertex = index_to_vertex.at(i);

    if( fiedlerVector[i] < median )
      partition[vertex] = true;
    else
      partition[vertex] = false;
  }

  std::vector<Simplex> simplices( K.begin(), K.end() );

  auto itLeft = std::stable_partition( simplices.begin(), simplices.end(),
    [&partition] ( const Simplex& s )
    {
      // All vertices of the simplex need to be part of the same
      // partition with respect to the matrix.
      return s.size() == IndexType( std::count_if( s.begin(), s.end(),
        [&partition] ( VertexType v )
        {
          return partition.at(v);
        }
      ) );
    }
  );

  auto itRight = std::stable_partition( itLeft, simplices.end(),
    [&partition] ( const Simplex& s )
    {
      // All vertices of the simplex need to be part of the same
      // partition with respect to the matrix.
      return s.size() == IndexType( std::count_if( s.begin(), s.end(),
        [&partition] ( VertexType v )
        {
          return !partition.at(v);
        }
      ) );
    }
  );

  std::vector<SimplicialComplex> complexes;
  complexes.push_back( SimplicialComplex( simplices.begin(), itLeft ) );
  complexes.push_back( SimplicialComplex( itLeft, itRight ) );

  return complexes;

#else
  (void) K;
  return {};
#endif
}

} // namespace topology

} // namespace aleph

#endif
