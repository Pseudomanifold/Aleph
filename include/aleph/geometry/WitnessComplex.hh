#ifndef ALEPH_GEOMETRY_WITNESS_COMPLEX_HH__
#define ALEPH_GEOMETRY_WITNESS_COMPLEX_HH__

#include <algorithm>
#include <limits>
#include <iterator>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

#include <aleph/geometry/RipsExpander.hh>

#include <aleph/geometry/distances/Traits.hh>

#include <aleph/math/SymmetricMatrix.hh>

#include <aleph/topology/Simplex.hh>
#include <aleph/topology/SimplicialComplex.hh>

#include <aleph/topology/filtrations/Data.hh>

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

  @param dimension Maximum dimension for expanding the witness complex
                   after obtaining its edges. The expansion process is
                   going to use the maximum possible dimension if this
                   parameter is not specified by the client.

  @param nu        \f$\nu\f$-parameter as defined in the paper. The
                   parameter controls which distance threshold is used
                   for creating an edge.\n

                   More precisely, an edge is created if the following
                   holds:\n
                   \n
                    \f[
                     \max(\mathrm{dist}_{a,i}\mathrm{dist}_{b,i}) \leq R + m_i
                    \f]
                   \n
                   Here, \f$m_i\f$ refers to the \f$\nu\f$th
                   smallest element of the distance matrix.

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
  unsigned dimension = 0,
  unsigned nu = 2,
  typename Distance::ResultType R = typename Distance::ResultType(),
  Distance /* distance */ = Distance() ) -> topology::SimplicialComplex< topology::Simplex<typename Distance::ResultType, typename std::iterator_traits<InputIterator>::value_type> >
{
  using IndexType         = typename std::iterator_traits<InputIterator>::value_type;
  using VertexType        = IndexType;
  using DataType          = typename Distance::ResultType;
  using Traits            = aleph::geometry::distances::Traits<Distance>;
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

  // Distance matrix between a set of $n$ landmarks (cols) and $N$ data
  // points (rows). Note that I transposed the matrix because accessing
  // the columns is faster that way (and will be required later on).
  std::vector< std::vector<DataType> > D;
  D.reserve( N );

  Distance dist;
  Traits traits;

  for( std::size_t j = 0; j < N; j++ )
  {
    std::vector<DataType> distances;
    distances.reserve( n );

    auto&& point = container[j];

    for( std::size_t i = 0; i < n; i++ )
    {
      auto&& landmark = container[ landmarkIndices.at(i) ];

      distances.emplace_back( traits.from( dist( landmark.begin(), point.begin(), d ) ) );
    }

    D.push_back( distances );
  }

  // Get smallest entries of the distance matrix. This is required for
  // deciding whether a specific edge is valid or not, with respect to
  // the given parameters.

  std::vector<DataType> smallest( N );

  if( nu != 0 )
  {
    for( std::size_t col = 0; col < N; col++ )
    {
      std::vector<DataType> column = D[col];

      std::nth_element( column.begin(), column.begin() + nu - 1, column.end() );
      smallest[col] = column.at( nu - 1 );
    }
  }

  // -------------------------------------------------------------------
  //
  // Records the appearance times of each potential edge in the witness
  // complex and creates the valid edges.

  std::vector<Simplex> simplices;

  aleph::math::SymmetricMatrix<DataType> M( n );

  for( std::size_t i = 0; i < n; i++ )
  {
    simplices.push_back( Simplex( static_cast<VertexType>(i) ) );

    for( std::size_t j = i+1; j < n; j++ )
    {
      auto min = std::numeric_limits<DataType>::max();

      for( std::size_t k = 0; k < N; k++ )
      {
        if( std::max( D[k][i], D[k][j] ) <= R + smallest.at(k) )
          min = std::min( min, std::max( D[k][i], D[k][j] ) );
      }

      if( min != std::numeric_limits<DataType>::max() )
      {
        auto u = static_cast<VertexType>(i);
        auto v = static_cast<VertexType>(j);

        simplices.push_back( Simplex( {u,v}, min ) );
      }
    }
  }

  aleph::geometry::RipsExpander<SimplicialComplex> ripsExpander;

  SimplicialComplex K = SimplicialComplex( simplices.begin(), simplices.end() );
  SimplicialComplex L = ripsExpander( K, dimension == 0 ? static_cast<unsigned>( d + 1 ) : dimension );
  L                   = ripsExpander.assignMaximumWeight( L );

  L.sort( aleph::topology::filtrations::Data<Simplex>() );
  return L;
}

/**
  Generates a random set of landmarks for use with the witness complex.
  Essentially, this function merely generates a set of *random indices*
  based on a random shuffle operation. The results of this are saved in
  an output iterator.

  @param n      Number of points
  @param k      Number of landmarks
  @param result Output iterator for storing the results
*/

template <class T, class OutputIterator> void generateRandomLandmarks( T n, T k, OutputIterator result )
{
  std::random_device rd;
  std::mt19937 rng( rd() );

  std::vector<T> indices( n );

  using DifferenceType = typename std::vector<T>::difference_type;

  std::iota( indices.begin(), indices.end(), T() );
  std::shuffle( indices.begin(), indices.end(), rng );
  std::copy( indices.begin(), indices.begin() + static_cast<DifferenceType>(k), result );
}

/**
  Generates a set of landmarks for the witness complex, using the
  max-min strategy. Given a distance measure, a new landmark will
  be chosen so as to *maximize* the *minimum distance* to the set
  of selected landmarks. An output iterator is used to report the
  indices of the selected landmarks.

  @param container Container that stores the input data
  @param n         Number of landmarks to select
  @param result    Output iterator for storing the results
  @param distance  Distance measure. This parameter may be specified
                   to permit template type deduction.
*/

template <
  class Distance,
  class Container,
  class OutputIterator
> void generateMaxMinLandmarks( const Container& container, std::size_t n, OutputIterator result, Distance distance = Distance() )
{
  if( n > container.size() )
    throw std::out_of_range( "Number of landmarks is out of range" );

  using SizeType = decltype( container.size() );

  std::random_device rd;
  std::mt19937 rng( rd() );

  std::uniform_int_distribution<SizeType> distribution( SizeType(0), container.size() - 1 );

  std::vector<SizeType> indices;
  indices.reserve( n );

  indices.emplace_back( distribution( rng ) );

  using DataType = typename Distance::ResultType;
  auto N         = container.size();
  auto d         = container.dimension();

  while( indices.size() < n )
  {
    auto index = SizeType(0);
    auto max   = std::numeric_limits<DataType>::lowest();

    for( SizeType i = 0; i < N; i++ )
    {
      auto min = std::numeric_limits<DataType>::max();

      for( auto&& landmarkIndex : indices )
      {
        auto dist = distance( container[i].begin(), container[landmarkIndex].begin(), d );
        min       = std::min( min, dist );
      }

      if( min > max )
      {
        max   = min;
        index = i;
      }
    }

    indices.push_back( index );
  }

  std::copy( indices.begin(), indices.end(), result );
}

} // namespace geometry

} // namespace aleph

#endif
