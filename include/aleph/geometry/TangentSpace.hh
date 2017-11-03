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

#ifdef ALEPH_WITH_EIGEN

class LocalTangentSpace
{
};

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
  };

  template <class Container> std::vector<T> operator()( const Container& container )
  {
    std::vector<double> curvature;
    curvature.reserve( container.size() );

    localTangentSpaces( container );

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
    auto k = 10;

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
      Matrix M = Matrix::Zero(k,d);

      // Centroid calculation ------------------------------------------

      Vector centroid  = Vector::Zero(1,d);

      for( std::size_t j = 0; j < indices[i].size(); j++ )
      {
        auto&& neighbourIndex = indices[i][j];
        auto p                = container[neighbourIndex];
        Vector v              = Vector::Zero(1,d);

        // copy (and transform!) the vector; there's an implicit type
        // conversion going on here
        for( std::size_t l = 0; l < d; l++ )
          v(l) = p[l];

        centroid += v;
        M.row(j)  = v;
      }

      centroid /= static_cast<T>( indices.size() );

      // Coordinate matrix setup ---------------------------------------

      M = M.rowwise() - centroid;

      Eigen::JacobiSVD<Matrix> svd( M, Eigen::ComputeThinV );

      LocalTangentSpace lts;
      lts.tangents = Matrix::Zero( d, d - 1 );

      // The singular vectors of all but the *smallest* singular value
      // form the tangential directions of the tangent space.

      auto&& V = svd.matrixV();

      for( std::size_t j = 0; j < d - 1; j++ )
        lts.tangents.col(j) = V.col(j);

      lts.normal = Matrix::Zero( 1, d );
      lts.normal = V.col(d-1);

      localTangentSpaces.push_back( lts );

      // TODO: calculate raw approximation (reconstruction) error by
      // assessing how well the space fits the original data set
    }

    return localTangentSpaces;
  }
};

#endif

} // namespace geometry

} // namespace aleph

#endif
