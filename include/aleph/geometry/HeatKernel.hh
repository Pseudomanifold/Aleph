#ifndef ALEPH_GEOMETRY_HEAT_KERNEL_HH__
#define ALEPH_GEOMETRY_HEAT_KERNEL_HH__

#include <aleph/config/Eigen.hh>

#ifdef ALEPH_WITH_EIGEN
  #include <Eigen/Core>
  #include <Eigen/Eigenvalues>
#endif

#include <aleph/math/KahanSummation.hh>

#include <algorithm>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <vector>

#include <cmath>

#define THROW_EIGEN_REQUIRED_ERROR()\
{\
  auto message =  std::string( __FILE__ )           \
                + std::string( ":" )                \
                + std::to_string( __LINE__ )        \
                + std::string( " in " )             \
                + std::string( __PRETTY_FUNCTION__ )\
                + std::string( ":" )                \
                + std::string( " Eigen is required for this function to work properly" );\
  \
  throw std::runtime_error( message );\
}

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

  Matrix W = Matrix::Zero( n, n );

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

/**
  Calculates the weighhted Laplacian matrix of a given simplicial
  complex and returns it.

  @param K Simplicial complex

  @returns Weighted Laplacian matrix. The indices of rows and columns
           follow the order of the vertices in the complex.
*/

template <class SimplicialComplex> auto weightedLaplacianMatrix( const SimplicialComplex& K ) -> Eigen::Matrix<typename SimplicialComplex::ValueType::DataType, Eigen::Dynamic, Eigen::Dynamic>
{
  auto W          = weightedAdjacencyMatrix( K );
  using Matrix    = decltype(W);
  using IndexType = typename Matrix::Index;

  Matrix L = Matrix::Zero( W.rows(), W.cols() );

  auto V = W.rowwise().sum();

  for( IndexType i = 0; i < V.size(); i++ )
    L(i,i) = V(i);

  return L - W;
}

#endif

/**
  @class HeatKernel
  @brief Calculates the heat kernel for simplicial complexes

  This class acts as a query functor for the heat kernel values of
  vertices in a weighted simplicial complex. It will pre-calculate
  the heat matrix and permit queries about the progression of heat
  values for *all* vertices for some time \f$t\f$.
*/

class HeatKernel
{
public:

  using T = double;

#ifdef ALEPH_WITH_EIGEN
  using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  using Vector = Eigen::Matrix<T, 1, Eigen::Dynamic>;

#if EIGEN_VERSION_AT_LEAST(3,3,0)
  using IndexType  = Eigen::Index;
#else
  using IndexType  = typename Matrix::Index;
#endif

#else
  // This declares a fallback index type in case Eigen is not available,
  // making sure that the interface of the class remains the same.
  using IndexType = unsigned;
#endif

  /**
    Constructs a heat kernel from a given simplicial complex. Afterwards,
    the functor will be ready for queries.

    @param K Simplicial complex
  */

  template <class SimplicialComplex> HeatKernel( const SimplicialComplex& K )
  {
#ifdef ALEPH_WITH_EIGEN

    auto L = weightedLaplacianMatrix( K );

    Eigen::SelfAdjointEigenSolver< decltype(L) > solver;
    solver.compute( L );

    auto&& eigenvalues  = solver.eigenvalues(). template cast<T>();
    auto&& eigenvectors = solver.eigenvectors().template cast<T>();

    _eigenvalues.reserve( std::size_t( eigenvalues.size() ) );
    _eigenvectors.reserve( std::size_t( eigenvectors.size() ) );

    using IndexType_ = typename decltype(L)::Index;

    // If configured, skip both the first eigenvector and the first
    // eigenvalue because they do not contribute anything later on.
    for( IndexType_ i = _skip ? 1 : 0; i < eigenvalues.size(); i++ )
      _eigenvalues.push_back( eigenvalues(i) );

    for( IndexType_ i = _skip ? 1 : 0; i < eigenvectors.cols(); i++ )
      _eigenvectors.push_back( eigenvectors.col(i) );

#else
  (void) K;
  (void) L;

  THROW_EIGEN_REQUIRED_ERROR();
#endif

  }

  /**
    Evaluates the heat kernel for *all* vertices at a given time \f$t\f$
    and returns the resulting values. This function is guaranteed to be
    more efficient than calling the per-element functions repeatedly.
  */

  std::vector<T> operator()( T t )
  {
#ifdef ALEPH_WITH_EIGEN

    Vector result = Vector();

    for( std::size_t k = 0; k < _eigenvalues.size(); k++ )
    {
      auto&& lk  = std::exp( -t * _eigenvalues[k] );
      auto&& uk = _eigenvectors[k];

      result += lk * uk * uk.transpose();
    }

    return std::vector<T>( result.data(), result.data() + result.size() );

#else
    (void) t;

    THROW_EIGEN_REQUIRED_ERROR();
#endif
  }

  /**
    Evaluates the heat kernel for two vertices \f$i\f$ and \f$j\f$ at
    a given time \f$t\f$ and returns the result.
  */

  T operator()( IndexType i, IndexType j, T t )
  {
#ifdef ALEPH_WITH_EIGEN

    aleph::math::KahanSummation<T> result = T();

    for( std::size_t k = 0; k < _eigenvalues.size(); k++ )
    {
      auto&& lk  = std::exp( -t * _eigenvalues[k] );
      auto&& uik = _eigenvectors[k](i);
      auto&& ujk = _eigenvectors[k](j);

      result += lk * uik * ujk;
    }

    return result;

#else
  (void) i;
  (void) j;
  (void) t;

  THROW_EIGEN_REQUIRED_ERROR();
#endif
  }

  /**
    Calculates the auto-diffusion for a given vertex \f$i\f$ and a given
    time \f$t\f$ and returns it.
  */

  T operator()( IndexType i, T t )
  {
#ifdef ALEPH_WITH_EIGEN

    // Note that this function could have been implemented in terms of
    // operator(i,j,t), but this implementation is a *little* bit more
    // efficient as it defines the multiplication explicitly.

    aleph::math::KahanSummation<T> result = T();

    for( std::size_t k = 0; k < _eigenvalues.size(); k++ )
    {
      auto&& lk  = std::exp( -t * _eigenvalues[k] );
      auto&& uik = _eigenvectors[k](i);

      result += lk * uik * uik;
    }

    return result;

#else
  (void) i;
  (void) t;

  THROW_EIGEN_REQUIRED_ERROR();
#endif
  }

  /**
    Calculates the *trace* of the heat kernel for a given time \f$t\f$
    and returns it.
  */

  T trace( T t ) const
  {
    aleph::math::KahanSummation<T> result = T();

    for( auto&& eigenvalue : _eigenvalues )
      result += std::exp( -t * eigenvalue );

    return result;
  }

  /**
    Calculates the *determinant* of the heat kernel for a given time
    \f$t\f$ and returns it.
  */

  T determinant( T t) const
  {
    T result = T();

    for( auto&& eigenvalue : _eigenvalues )
      result = result * std::exp( -t * eigenvalue );

    return result;
  }

  // Sampling intervals ------------------------------------------------

  /**
    Uses a heuristic to determine a sampling interval for the time
    parameter \f$t\f$ of the heat kernel. This heuristic was first
    described by Sun et al. in their paper *A Concise and Provably
    Informative Multi-Scale Signature Based on Heat Diffusion*.

    @param n Number of sampling points
    @returns Vector of sampling points
  */

  std::vector<T> logarithmicSamplingInterval( unsigned n ) const
  {
    auto t_min  = 4 * std::log( 10 ) / _eigenvalues.back();
    auto t_max  = 4 * std::log( 10 ) / ( _skip ? _eigenvalues.front() : *( _eigenvalues.begin() + 1 ) );
    auto offset = ( std::log( t_max ) - std::log( t_min ) ) / ( n - 1 );

    std::vector<T> samples;
    samples.reserve( n );

    for( unsigned i = 0; i < n; i++ )
      samples.push_back( std::log( t_min ) + i * offset );

    std::transform( samples.begin(), samples.end(),
                    samples.begin(),
                    [] ( const T x )
                    {
                      return std::pow( std::exp(1), x );
                    } );

    return samples;
  }

  // Configuration -----------------------------------------------------

  void setSkip( bool value = true ) { _skip = value; }
  bool skip() const noexcept        { return _skip;  }

private:

  /** If set, skips the first eigenvector and eigenvalue */
  bool _skip = false;

  /**
    Stores the eigenvalues of the heat matrix, or, more precisely, the
    eigenvalues of the Laplacian. They will be used for the evaluation
    of the heat kernel.
  */

  std::vector<T> _eigenvalues;

#ifdef ALEPH_WITH_EIGEN

  /** Stores the eigenvectors of the heat matrix */
  std::vector<Vector> _eigenvectors;

  /**
    Heat matrix; will be created automatically upon constructing this
    functor class.
  */

  Matrix _H;
#endif

};

} // namespace geometry

} // namespace aleph

#endif
