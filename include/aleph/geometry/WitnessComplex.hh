#ifndef ALEPH_GEOMETRY_WITNESS_COMPLEX_HH__
#define ALEPH_GEOMETRY_WITNESS_COMPLEX_HH__

#include <algorithm>
#include <limits>
#include <iterator>
#include <vector>

#include <aleph/geometry/distances/Traits.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

namespace aleph
{

namespace geometry
{

template <class Distance, class Container, class InputIterator> auto buildWitnessComplex(
  const Container& container,
  InputIterator begin, InputIterator end,
  unsigned nu = 2,
  typename Distance::ResultType R = typename Distance::ResultType() ) -> topology::SimplicialComplex< topology::Simplex<typename Distance::ResultType, unsigned> >
{
  using IndexType         = typename std::iterator_traits<InputIterator>::value_type;
  using DataType          = typename Distance::ResultType;
  using Traits            = aleph::distances::Traits<Distance>;
  using VertexType        = unsigned; // TODO: make configurable
  using Simplex           = topology::Simplex<DataType, VertexType>;
  using SimplicialComplex = topology::SimplicialComplex<Simplex>;

  // These are only the *indices* of the landmarks, with respect to the
  // underlying point cloud.
  std::vector<IndexType> landmarkIndices( begin, end );

  auto n = landmarkIndices.size();
  auto N = container.size();
  auto d = container.dimension();

  // Much of the behaviour below will be undefined if we permit such
  // situations to occur.
  if( n == 0 || N == 0 )
    return {};

  // Distance matrix between a set of $n$ landmarks (rows) and $N$ data
  // points.
  std::vector< std::vector<DataType> > D;
  D.reserve( n );

  Distance dist;
  Traits traits;

  for( std::size_t i = 0; i < n; i++ )
  {
    std::vector<DataType> distances;
    distances.reserve( N );

    auto&& landmark = container[ landmarkIndices.at(i) ];

    for( std::size_t j = 0; j < N; j++ )
    {
      auto&& point = container[j];

      distances.emplace_back( traits.from( dist( landmark.begin(), point.begin(), d ) ) );
    }

    D.push_back( distances );
  }

  // Records the appearance times of each potential edge in the witness
  // complex.
  //
  // TODO: this should become a symmetric matrix
  std::vector< std::vector<DataType> > M( n, std::vector<DataType>( n ) );

  for( std::size_t i = 0; i < n; i++ )
  {
    for( std::size_t j = i+1; j < n; j++ )
    {
      auto min = std::numeric_limits<DataType>::max();

      for( std::size_t k = 0; k < N; k++ )
        min = std::min( min, std::max( D[i][k], D[j][k] ) );

      M[i][j] = min;
      M[j][i] = min;
    }
  }

  // Get smallest entries of the distance matrix. This is required for
  // deciding whether a specific edge is valid or not, with respect to
  // the given parameters.

  std::vector<DataType> smallest( N );

  if( nu != 0 )
  {
    for( std::size_t col = 0; col < N; col++ )
    {
      // FIXME: getting the column like this is extremely wasteful;
      // would it not be nicer to store the values differently?
      std::vector<DataType> column( n );
      for( std::size_t i = 0; i < n; i++ )
        column[i] = D[i][col];

      std::nth_element( column.begin(), column.begin() + nu - 1, column.end() );
      smallest[col] = column.at( nu - 1 );
    }
  }

  auto max = *std::max_element( smallest.begin(), smallest.end() );

  std::vector<Simplex> simplices;

  for( std::size_t i = 0; i < n; i++ )
  {
    simplices.push_back( Simplex( static_cast<VertexType>(i) ) );

    for( std::size_t j = i+1; j < n; j++ )
    {
      // Skip pairs that cannot possibly give rise to an edge because of
      // their distance to each other.
      if( M[i][j] > R + max )
        continue;

      for( std::size_t col = 0; col < N; col++ )
      {
        if( M[i][j] <= R + smallest.at(col) )
        {
          auto u = static_cast<VertexType>(i);
          auto v = static_cast<VertexType>(j);

          simplices.push_back( Simplex( {u,v}, M[i][j] ) );
          break;
        }
      }
    }
  }

  return SimplicialComplex( simplices.begin(), simplices.end() );
}

} // namespace geometry

} // namespace aleph

#endif
