#ifndef ALEPH_GEOMETRY_TANGENT_SPACE_HH__
#define ALEPH_GEOMETRY_TANGENT_SPACE_HH__

#include <aleph/config/Eigen.hh>

#include <aleph/geometry/BruteForce.hh>
#include <aleph/geometry/FLANN.hh>

#include <aleph/geometry/distances/Euclidean.hh>

#ifdef ALEPH_WITH_EIGEN
  #include <Eigen/Core>
  #include <Eigen/SVD>
#endif

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

class TangentSpace
{
public:

  using T        = double;

  using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  using Vector = Eigen::Matrix<T, 1, Eigen::Dynamic>;

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

    std::vector<std::size_t> indices;
  };

  template <class Container> std::vector<T> operator()( const Container& container )
  {
    std::vector<double> curvature;
    curvature.reserve( container.size() );

    auto lts = localTangentSpaces( container );
    fitSpheres( container, lts );

    return curvature;
  }

  template <class Container> std::vector<LocalTangentSpace> localTangentSpaces( const Container& container )
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

    // FIXME: make configurable
    unsigned k = 10;

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
        auto p                = container[neighbourIndex];
        Vector v              = Vector::Zero(1, Index(d) );

        // copy (and transform!) the vector; there's an implicit type
        // conversion going on here
        for( std::size_t l = 0; l < d; l++ )
          v( Index(l) ) = p[l];

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

      lts.normal   = Matrix::Zero( Index(1), Index(d) );
      lts.normal   = V.col( Index(d-1) );
      lts.position = getPosition( container, i );
      lts.indices  = indices[i];

      localTangentSpaces.push_back( lts );

      // TODO: calculate raw approximation (reconstruction) error by
      // assessing how well the space fits the original data set
    }

    return localTangentSpaces;
  }

  template <class Container> void fitSpheres( const Container& container,
                                              const std::vector<LocalTangentSpace>& localTangentSpaces )
  {
    using namespace detail;

    for( auto&& lts : localTangentSpaces )
    {
      auto d = Index( container.dimension() );
      auto w = 1.0;

      Matrix A = Matrix::Zero( d+2, d+2 );
      Vector b = Vector::Zero( 1,   d+2 );

      auto&& indices = lts.indices;

      for( auto&& index : indices )
      {
        auto neighbour           = getPosition( container, index );
        auto squaredNeigbourNorm = neighbour.squaredNorm();
        auto localFeatureSize    = 1.0; // FIXME: how to obtain estimates?
        auto _beta               = 1.0; // FIXME: how to obtain estimates?
        w                        = phi( ( lts.position - neighbour ).norm() / localFeatureSize );

        A(   0,   0) += w;
        A( d+1,   0) += w * squaredNeigbourNorm;
        A( d+1, d+1) += w * squaredNeigbourNorm * squaredNeigbourNorm;

        for( Index(i) = 1; i < d+1; i++ )
        {
          A(  i,   i) += w * ( neighbour(i-1)*neighbour(i-1) + 1 ) * _beta;
          A(  i,   0) += w * ( neighbour(i-1) );
          A(d+1,   i) += w * ( neighbour(i-1)*squaredNeigbourNorm + 2 * _beta * neighbour(i-1) );
          A(d+1, d+1) += w * ( 4*neighbour(i-1)*neighbour(i-1) ) * _beta;

          // re-establish symmetry
          A(  0,   i)  = A(i,0);
          A(  i, d+1)  = A(d+1, i);

          b(i  ) +=       _beta * w * lts.normal(i-1);
          b(d+1) += 2.0 * _beta * w * lts.normal(i-1) * neighbour(i-1);

          for( Index(j) = i+1; j < d+1; j++ )
          {
            A(j, i) += w * neighbour(i-1) * neighbour(j-1);
            A(i, j)  = A(j,i);
          }
        }
      }

      std::cout << A << "," << b << "\n";
    }
  }

private:

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

} // namespace geometry

} // namespace aleph

#endif
