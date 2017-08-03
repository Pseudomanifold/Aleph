#ifndef ALEPH_GEOMETRY_WITNESS_COMPLEX_HH__
#define ALEPH_GEOMETRY_WITNESS_COMPLEX_HH__

#include <algorithm>
#include <limits>
#include <iterator>
#include <vector>

#include <aleph/geometry/RipsExpander.hh>

#include <aleph/geometry/distances/Traits.hh>

#include <aleph/math/SymmetricMatrix.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

namespace aleph
{

namespace geometry
{

/**
  Builds a witness complex from a given container. This requires a set
  of *landmarks*. Other configuration options influence how a new edge
  will be created from the data.

  If you call this function with its barest minimum parameters by only
  specifying a container and a set of landmarks, the resulting complex
  will automatically adjust to your data and follows the definition of
  the paper:

    > Topological estimation using witness complexes\n
    > Vin de Silva and Gunnar Carlsson\n
    > Eurographics Symposium on Point-Based Graphics, 2004

  The other parameters, in particular \p R, permit tuning the results,
  thereby given the complex more "slack" when creating edges. However,
  this also increases the size of the complex.

  @param container Container for which to calculate the witness complex

  @param begin     Input iterator to begin of landmark range; landmarks
                   must be specified as *indices*. Each index has to be
                   valid for \p container.

  @param end       Input iterator to end of landmark range

  @param nu        \f$\nu\f$-parameter as specified in the paper. This
                   parameter controls which distance threshold is used
                   for creating an edge.\n

                   More precisely, an edge is created if the following
                   holds:\n

                   \f[
                    \max(\mathrm{dist}_{a,i}\mathrm{dist}_{b,i}) \leq R + m_i
                   \f]
                   Here, \f$m_i\f$ refers to the \f$\nu\f$-th smallest
                   element of the distance matrix.

  @param R         Maximum radius parameter for determining whether to
                   create an edge or not. The default value makes sure
                   that *only* the distances between the landmarks and
                   the points will be used for edge creation.

  @param distance  Distance functor, e.g. the Euclidean distance, that
                   is used for calculating distances between landmarks
                   and data points. This parameter is provided to make
                   it possible for the compiler to detect the template
                   parameter \p distance automatically.

  @returns         Witness complex of the given container. Notice that
                   the complex is stored as a simplicial complex whose
                   data type and index type are derived from the input
                   data. The data of a simplex contains the *smallest*
                   threshold for which they appear in the complex. The
                   complex will be sorted according to this value.
*/


template <
  class Distance,
  class Container,
  class InputIterator
> auto buildWitnessComplex(
  const Container& container,
  InputIterator begin,
  InputIterator end,
  unsigned nu = 2,
  typename Distance::ResultType R = typename Distance::ResultType(),
  Distance /* distance */ = Distance() ) -> topology::SimplicialComplex< topology::Simplex<typename Distance::ResultType, typename std::iterator_traits<InputIterator>::value_type> >
{
  using IndexType         = typename std::iterator_traits<InputIterator>::value_type;
  using VertexType        = IndexType;
  using DataType          = typename Distance::ResultType;
  using Traits            = aleph::distances::Traits<Distance>;
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

  aleph::math::SymmetricMatrix<DataType> M( n );

  for( std::size_t i = 0; i < n; i++ )
  {
    for( std::size_t j = i+1; j < n; j++ )
    {
      auto min = std::numeric_limits<DataType>::max();

      for( std::size_t k = 0; k < N; k++ )
        min = std::min( min, std::max( D[i][k], D[j][k] ) );

      M(i,j) = min;
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
      if( M(i,j) > R + max )
        continue;

      for( std::size_t col = 0; col < N; col++ )
      {
        if( M(i,j) <= R + smallest.at(col) )
        {
          auto u = static_cast<VertexType>(i);
          auto v = static_cast<VertexType>(j);

          // FIXME: this is incorrect; there could be edges with
          // a smaller weight
          simplices.push_back( Simplex( {u,v}, M(i,j) ) );
          break;
        }
      }
    }
  }

  aleph::geometry::RipsExpander<SimplicialComplex> ripsExpander;

  SimplicialComplex K = SimplicialComplex( simplices.begin(), simplices.end() );
  SimplicialComplex L = ripsExpander( K, static_cast<unsigned>(n) ); // TODO: make dimension configurable?
  L                   = ripsExpander.assignMaximumWeight( L );

  return L;
}

} // namespace geometry

} // namespace aleph

#endif
