#ifndef ALEPH_GEOMETRY_TANGENT_SPACE_HH__
#define ALEPH_GEOMETRY_TANGENT_SPACE_HH__

#include <aleph/config/Eigen.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/FLANN.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#include <aleph/math/AlgebraicSphere.hh>
#include <aleph/math/KahanSummation.hh>

#ifdef ALEPH_WITH_EIGEN
  #include <Eigen/Core>
  #include <Eigen/Cholesky>
  #include <Eigen/SVD>
#endif

#include <cassert>

#include <vector>

namespace aleph
{

namespace geometry
{

namespace detail
{

/**
  Model of a smooth decreasing weight function according to the
  original paper *Algebraic Point Set Surfaces* by Guennebaud &
  Gross.
*/

template <class T> T phi( T x )
{
  return x < 1 ? std::pow( 1 - x*x, T(4) ) : T();
}

} // namespace detail

#ifdef ALEPH_WITH_EIGEN

// Previous versions of Eigen have a bug that occurs when mixing dynamic
// and fixed-sized vectors:
//
//   http://eigen.tuxfamily.org/bz/show_bug.cgi?id=654
//
// Until a workaround has been identified, tangent space estimation will
// not be enabled for older versions.
#if EIGEN_VERSION_AT_LEAST(3,3,0)

class TangentSpace
{
public:
  using T        = double;

  using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  using Vector = Eigen::Matrix<T, 1, Eigen::Dynamic>;

  using Sphere = math::AlgebraicSphere<T>;

#if EIGEN_VERSION_AT_LEAST(3,3,0)
  using Index = Eigen::Index;
#else
  using Index = typename Matrix::Index;
#endif

  struct LocalTangentSpace
  {
    Matrix tangents;
    Vector normal;
    Vector position;

    T localFeatureSize;
    std::vector<std::size_t> indices;
  };

  template <class Container> std::vector<T> operator()( const Container& container, unsigned k )
  {
    std::vector<double> curvature;
    curvature.reserve( container.size() );

    auto lts     = localTangentSpaces( container, k );
    auto spheres = fitSpheres( container, lts );

    std::transform( spheres.begin(), spheres.end(), std::back_inserter( curvature ),
      [] ( const Sphere& sphere )
      {
        return sphere.meanCurvature();
      }
    );

    return curvature;
  }

private:

  /**
    Given a container and a number for determining the local
    neighbourhood of points, this function estimates (!) the
    tangent space structure around every point, resulting in
    a set of normal vectors and tangent vectors.

    @param container Container
    @param k         Local neighbourhood size
  */

  template <class Container> std::vector<LocalTangentSpace> localTangentSpaces( const Container& container, unsigned k )
  {
    using ElementType = typename Container::ElementType;
    using Distance    = distances::Euclidean<ElementType>;

#ifdef ALEPH_WITH_FLANN
    using NearestNeighbours = FLANN<Container, Distance>;
#else
    using NearestNeighbours = BruteForce<Container, Distance>;
#endif

    NearestNeighbours nearestNeighbours( container );
    using IndexType = typename NearestNeighbours::IndexType;

    std::vector< std::vector<IndexType> > indices;
    std::vector< std::vector<ElementType> > distances;

    nearestNeighbours.neighbourSearch( k,
                                       indices,
                                       distances );

    std::vector<LocalTangentSpace> localTangentSpaces;

    auto n = container.size();
    auto d = container.dimension();

    for( std::size_t i = 0; i < n; i++ )
    {
      // This coordinate matrix will contain the differences to the
      // centroid coordinate. The matrix will be transformed via an
      // SVD.
      Matrix M = Matrix::Zero( Index(k), Index(d) );

      // Centroid calculation ------------------------------------------

      Vector centroid  = Vector::Zero(1, Index(d) );

      for( std::size_t j = 0; j < indices[i].size(); j++ )
      {
        auto&& neighbourIndex = indices[i][j];
        auto v                = getPosition( container, neighbourIndex );

        centroid          += v;
        M.row( Index(j) )  = v;
      }

      centroid /= static_cast<T>( indices.size() );

      // Coordinate matrix setup ---------------------------------------

      M = M.rowwise() - centroid;

      Eigen::JacobiSVD<Matrix> svd( M, Eigen::ComputeThinV );

      LocalTangentSpace lts;
      lts.tangents = Matrix::Zero( Index(d), Index(d - 1) );

      // The singular vectors of all but the *smallest* singular value
      // form the tangential directions of the tangent space.

      auto&& V = svd.matrixV();

      for( Index j = 0; j < Index( d - 1 ); j++ )
        lts.tangents.col(j) = V.col(j);

      lts.normal           = Matrix::Zero( Index(1), Index(d) );
      lts.normal           = V.col( Index(d-1) ).normalized();
      lts.position         = getPosition( container, i );
      lts.indices          = indices[i];

      // Take the *maximum distance* in which we can find all of the
      // neighbours as a *rough*  approximation to the local feature
      // size.
      lts.localFeatureSize
        = distances[i].empty() == false ?
            *std::max_element( distances[i].begin(), distances[i].end() )
          : T();

      localTangentSpaces.push_back( lts );

      // TODO:
      //
      // 1. Calculate raw approximation (reconstruction) error by
      //    assessing how well the space fits the original data
      //
      // 2. Make normal orientation consistent. I am unsure as to
      //    whether this will improve the results or not.
    }

    return localTangentSpaces;
  }

  template <class Container>
    std::vector<Sphere> fitSpheres( const Container& container,
                                    const std::vector<LocalTangentSpace>& localTangentSpaces )
  {
    using namespace detail;

    std::vector<Sphere> spheres;
    spheres.reserve( container.size() );

    for( auto&& lts : localTangentSpaces )
    {
      auto d   = Index( container.dimension() );
      Matrix A = Matrix::Zero( d+2, d+2 );
      Vector b = Vector::Zero( 1,   d+2 );

      auto&& indices = lts.indices;

      // Pre-processing --------------------------------------------------
      //
      // Choose a value for the beta parameter, based on the weighted
      // neighbourhood sizes of *all* points. This requires iterating
      // over all points prior to calculating anything else.

      T beta = T();

      {
        std::vector<T> W;
        std::vector<T> H;
        W.reserve( container.size() );
        H.reserve( container.size() );

        for( auto&& index : indices )
        {
          auto&& neighbour = getPosition( container, index );
          W.emplace_back( phi( ( lts.position - neighbour ).norm() / lts.localFeatureSize ) );
          H.emplace_back( lts.localFeatureSize );
        }

        // Sum of weights *before* applying the local feature size
        // multiplier; we need to save this result because it will
        // be required below.
        auto ws = math::accumulate_kahan_sorted( W.begin(), W.end(), T() );

        // Apply weights to local feature size estimates; afterwards we
        // can finally obtain the scaling factor from this weighted sum
        for( std::size_t i = 0; i < W.size(); i++ )
          W[i] = W[i] * H[i];

        // TODO: make initial guess for beta (10e6) configurable?
        auto h = math::accumulate_kahan_sorted( W.begin(), W.end(), T() ) / ws;
        beta   = 10e6 * h * h;
      }

      Index k  = Index( indices.size() );
      Matrix W = Matrix::Zero( (d+1)*k, (d+1)*k );
      Matrix D = Matrix::Zero( (d+1)*k, (d+2)   );
      Matrix c = Matrix::Zero( (d+1)*k, 1       );

      {
        Index i = Index();

        for( auto&& index : indices )
        {
          auto neighbour       = getPosition( container, index );
          auto w               = phi( ( lts.position - neighbour ).norm() / lts.localFeatureSize );
          W( i*(d+1),i*(d+1) ) = w;

          for( Index(j) = 0; j < d; j++ )
          {
            assert( W( i*(d+1)+j+1, i*(d+1)+j+1 ) == 0 );
            W( i*(d+1)+j+1, i*(d+1)+j+1 ) = beta * w;
          }

          ++i;
        }
      }

      {
        Index i = Index();

        for( auto&& index : indices )
        {
          auto neighbour       = getPosition( container, index );
          D( i*(d+1), 0 )      = 1.0;

          for( Index(j) = 0; j < d; j++ )
            D( i*(d+1), j+1 ) = neighbour(j);

          D( i*(d+1), d+1) = neighbour * neighbour.transpose();

          for( Index(j) = 0; j < d; j++ )
          {
            D( i*(d+1)+j+1, j+1 ) = 1;
            D( i*(d+1)+j+1, d+1 ) = 2 * neighbour(j);
            c( i*(d+1)+j+1, 0   ) = localTangentSpaces.at(index).normal(j);
          }

          ++i;
        }

      }

      for( auto&& index : indices )
      {
        auto neighbour           = getPosition( container, index );
        auto squaredNeigbourNorm = neighbour.squaredNorm();
        auto w                   = phi( ( lts.position - neighbour ).norm() / lts.localFeatureSize );

        A(   0,   0) += w;
        A( d+1,   0) += w * squaredNeigbourNorm;
        A( d+1, d+1) += w * squaredNeigbourNorm * squaredNeigbourNorm;

        for( Index(i) = 1; i < d+1; i++ )
        {
          A(  i,   i) += w * ( neighbour(i-1)*neighbour(i-1) + 1 ) * beta;
          A(  i,   0) += w * ( neighbour(i-1) );
          A(d+1,   i) += w * ( neighbour(i-1)*squaredNeigbourNorm + 2 * beta * neighbour(i-1) );
          A(d+1, d+1) += w * ( 4*neighbour(i-1)*neighbour(i-1) ) * beta;

          // re-establish symmetry
          A(  0,   i)  = A(i,0);
          A(  i, d+1)  = A(d+1, i);

          b(i  ) +=       beta * w * localTangentSpaces.at(index).normal(i-1);
          b(d+1) += 2.0 * beta * w * localTangentSpaces.at(index).normal(i-1) * neighbour(i-1);

          for( Index(j) = i+1; j < d+1; j++ )
          {
            A(j, i) += w * neighbour(i-1) * neighbour(j-1);
            A(i, j)  = A(j,i);
          }
        }

        // re-establish symmetry
        A(0, d+1) = A(d+1, 0);
      }

      // Solve the linear system ---------------------------------------
      //
      // The solution of the system Ax = b is used to obtain the
      // coefficients of the algebraic sphere.

      using Solver = Eigen::LDLT<Matrix>;
      Solver solver(A);

      Vector u = solver.solve( b.transpose() );

      spheres.emplace_back( Sphere( u.data(), u.data() + u.size() ) );
    }

    return spheres;
  }

  /**
    Auxiliary function for extracting and converting a position from
    a given container, storing it as a (mathematical) vector.
  */

  template <class Container> Vector getPosition( const Container& container, std::size_t i )
  {
    auto d   = container.dimension();
    auto p   = container[i];
    Vector v = Vector::Zero(1, Index(d) );

    // copy (and transform!) the vector; there's an implicit type
    // conversion going on here
    for( std::size_t l = 0; l < d; l++ )
      v( Index(l) ) = p[l];

    return v;
  }
};

#endif

#endif

} // namespace geometry

} // namespace aleph

#endif
