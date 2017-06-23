#ifndef ALEPH_CONTAINERS_POINT_CLOUD_HH__
#define ALEPH_CONTAINERS_POINT_CLOUD_HH__

#include <algorithm>
#include <fstream>
#include <initializer_list>
#include <iterator>
#include <stdexcept>
#include <string>

#include <cstddef>

#include <aleph/utilities/String.hh>

namespace aleph
{

namespace containers
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
    : PointCloud()
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

  friend void swap( PointCloud& pc1, PointCloud& pc2 ) noexcept
  {
    using std::swap;

    swap( pc1._points, pc2._points );
    swap( pc1._n,      pc2._n );
    swap( pc1._d,      pc2._d );
  }

  // Equality comparison -----------------------------------------------

  bool operator==( const PointCloud<T>& other ) const noexcept
  {
    return    _n == other._n
           && _d == other._d
           && std::equal( _points, _points + _n * _d, other._points );
  }

  // Attributes --------------------------------------------------------

  std::size_t size() const noexcept
  {
    return _n;
  }

  std::size_t dimension() const noexcept
  {
    return _d;
  }

  bool empty() const noexcept
  {
    return _n == 0;
  }

  // Point access ------------------------------------------------------

  // This is slightly evil. The function is not really "bit-wise"
  // constant because it still permits modifying the pointer that
  // is being returned.
  T* data() const noexcept
  {
    return _points;
  }

  /**
    Sets $i$th point of point cloud. Throws if the number of dimensions
    does not match the number of dimensions in the point cloud.
  */

  template <class InputIterator> void set( std::size_t i,
                                           InputIterator begin, InputIterator end )
  {
    if( i >= this->size() )
      throw std::runtime_error( "Invalid index" );

    auto distance = std::distance( begin, end );

    if( static_cast<std::size_t>( distance ) != this->dimension() )
      throw std::runtime_error( "Incorrect number of dimensions" );

    std::copy( begin, end, _points + this->dimension() * i );
  }

  /** @overload set() */
  void set( std::size_t i, const std::initializer_list<T>& il )
  {
    this->set( i, il.begin(), il.end() );
  }

  /**
    Gets the $i$th point of the point cloud. It will be stored via an
    output iterator. Incorrect indices will result in an exception.
  */

  template <class OutputIterator> void get( std::size_t i,
                                            OutputIterator result ) const
  {
    if( i >= this->size() )
      throw std::runtime_error( "Invalid index" );

    auto offset = i * this->dimension();

    std::copy( _points + offset, _points + offset + this->dimension(),
               result );
  }

  /**
    Stores the $i$th point of the point cloud in a vector and returns
    said vector. This interface is not very efficient but is provided
    as a fallback/simplified variant.
  */

  std::vector<T> operator[]( std::size_t i ) const
  {
    std::vector<T> p;

    this->get( i,
               std::back_inserter( p ) );

    return p;
  }

private:
  std::size_t _n; ///< Number of points
  std::size_t _d; ///< Dimension

  T* _points;
};

/**
  Loads a new point cloud from a file. The file is supposed to be in
  ASCII format. Each row must specify one item of the data set.  The
  different attributes of each item are assumed to be separated by a
  comma or white-space characters.
*/
template<class T> PointCloud<T> load( const std::string& filename )
{
  std::ifstream in( filename );

  if( !in )
    return PointCloud<T>();

  auto lines = std::count( std::istreambuf_iterator<char>( in ),
                           std::istreambuf_iterator<char>(),
                           '\n' );

  in.clear();
  in.seekg( 0 );

  if( lines <= 0 )
    return PointCloud<T>();

  std::size_t i = 0;
  std::size_t d = 0;
  std::size_t n = static_cast<std::size_t>( lines );

  std::string line;

  PointCloud<T> pointCloud;

  while( std::getline( in, line ) )
  {
    auto tokens = utilities::split( line, std::string( "[:;,[:space:]]+" ) );

    if( d == 0 )
    {
      d          = tokens.size();
      pointCloud = PointCloud<T>( n, d );
    }

    std::vector<T> coordinates;
    coordinates.reserve( d );

    for( auto&& token : tokens )
    {
      T coordinate = utilities::convert<T>( token );
      coordinates.push_back( coordinate );
    }

    pointCloud.set( i,
                    coordinates.begin(), coordinates.end() );

    ++i;
  }

  return pointCloud;
}

} // namespace containers

} // namespace aleph

#endif
