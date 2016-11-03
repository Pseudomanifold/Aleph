#ifndef ALEPH_CONTAINERS_POINT_CLOUD_HH__
#define ALEPH_CONTAINERS_POINT_CLOUD_HH__

#include <algorithm>

#include <cstddef>

namespace aleph
{

template <class T> class PointCloud
{
public:

  // Exporting the type of elements stored in the point cloud. This is
  // used by other algorithms to prevent casting.
  using ElementType = T;

  PointCloud()
    : _n( 0 )
    , _d( 0 )
    , _points( nullptr )
  {
  }

  PointCloud( std::size_t n, std::size_t d )
    : _n( n )
    , _d( d )
    , _points( new T[ _n * _d ] )
  {
    // Zero-initialization. This may not be the most efficient way of
    // doing it, in particular if a client has some data to pass, but
    // it ensures consistency.
    std::fill( _points, _points + _n * _d, T() );
  }

  PointCloud( const PointCloud& other )
    : _n( other._n )
    , _d( other._d )
    , _points( new T[ _n * _d ] )
  {
    std::copy( other._points, other._points + _n * _d, _points );
  }

  PointCloud( PointCloud&& other )
  {
    swap( *this, other );
  }

  PointCloud& operator=( PointCloud other )
  {
    swap( *this, other );
    return *this;
  }

  ~PointCloud()
  {
    delete[] _points;
  }

  friend void swap( PointCloud& pc1, PointCloud& pc2 )
  {
    using std::swap;

    swap( pc1._points, pc2._points );
    swap( pc1._n,      pc2._n );
    swap( pc1._d,      pc2._d );
  }

private:
  std::size_t _n; ///< Number of points
  std::size_t _d; ///< Dimension

  T* _points;
};

}

#endif
